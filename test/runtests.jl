using Test
using IgBench
import IgSim
using Random

const FIX = joinpath(@__DIR__, "fixtures")
const VFA = joinpath(FIX, "V.fasta")
const DFA = joinpath(FIX, "D.fasta")
const JFA = joinpath(FIX, "J.fasta")

@testset "matchers" begin
    @test allele_call_match("IGHV1-2*02", "IGHV1-2*02")
    @test allele_call_match("IGHV1-2*02", "IGHV1-2*02/IGHV1-2*01")
    @test allele_call_match("IGHV1-2*02_S1234", "IGHV1-2*02")
    @test gene_call_match("IGHV1-2*02", "IGHV1-2*01")
    @test !exact_call_match("IGHV1-2*02", "IGHV1-2*01")
    @test span_iou(Span(1, 10), Span(5, 15)) ≈ 6 / 15
    @test isnan(span_iou(EMPTY_SPAN, Span(1, 5)))
end

@testset "airr io roundtrip" begin
    rows = [
        CallRecord("s1", "ACGT", "V1*01", "D1*01", "J1*01", Span(1, 2), Span(3, 3), Span(4, 4)),
        CallRecord("s2", "TTTT", "V2*01", "", "J2*01"),
    ]
    path = joinpath(tempdir(), "igbench_airr_test.tsv.gz")
    write_airr_calls(path, rows)
    back = read_airr_calls(path)
    @test length(back) == 2
    @test back[1].v_call == "V1*01"
    @test back[1].v_span == Span(1, 2)
    @test back[2].d_call == ""
end

@testset "FakeAnnotator + metrics" begin
    gold = [
        CallRecord("a", "AAA", "V1*01", "D1*01", "J1*01", Span(1, 1), Span(2, 2), Span(3, 3)),
        CallRecord("b", "BBB", "V2*01", "", "J2*01", Span(1, 1), EMPTY_SPAN, Span(2, 2)),
    ]
    pred = [
        CallRecord("a", "AAA", "V1*01", "D1*01", "J1*01", Span(1, 1), Span(2, 2), Span(3, 3)),
        CallRecord("b", "BBB", "V2*01", "", "J9*01", Span(1, 1), EMPTY_SPAN, Span(2, 2)),
    ]
    m = evaluate(ExactCallAccuracy(), pred, gold)
    @test m.v == 1.0
    @test m.j == 0.5
    @test m.d == 1.0
    @test m.n == 2
    @test m.d_n == 1
end

@testset "SimGoldPanel + run_suite Fake" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    db = IgSim.load_germline(; v = VFA, d = DFA, j = JFA)
    hold = String[db.v[1].name]
    src = SimSource(; id = "toy", germline = gp, species = "toy_species",
                    holdout_v = hold, n = 8, seed = 1)
    panel = SimGoldPanel(src, :full; n = 8)
    data = load_panel(panel)
    @test length(data.sequences) == 8
    @test !isnothing(data.gold)

    tool = FakeAnnotator("perfect"; gold = data.gold)
    suite = BenchSuite("smoke";
                       panels = [panel],
                       tools = [tool],
                       compares = [CompareSpec("perfect", ":gold")],
                       timing = TimingSpec(warmup = 0, repeats = 1))
    outdir = mktempdir()
    store = DirectoryRunStore(outdir)
    result = run_suite(suite; mode = FullReportMode(timing_warmup = 0, timing_repeats = 1),
                       store, step = 42, tags = (; epoch = 1))
    @test result.step == 42
    @test any(r -> r["metric"] == "exact" && r["v"] == 1.0, result.metrics)
    @test isfile(joinpath(outdir, "metrics.jsonl"))
    @test isfile(joinpath(outdir, "metrics.json"))
    @test isfile(joinpath(outdir, "summary.md"))
    @test isfile(joinpath(outdir, "manifest.json"))
end

@testset "suite_from_manifest + diagnostic" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    db = IgSim.load_germline(; v = VFA, d = DFA, j = JFA)
    hold = String[db.v[1].name]
    # mini AIRR from sim gold
    src = SimSource(; id = "toy", germline = gp, species = "sp",
                    holdout_v = hold, n = 16, seed = 2)
    data = load_panel(SimGoldPanel(src, :full; n = 16))
    airr_path = joinpath(tempdir(), "igbench_mini.airr.tsv")
    write_airr_calls(airr_path, data.gold)
    man = DatasetManifest("mini";
                          sim = [src],
                          airr = [AirrSource(; id = "real_toy", path = airr_path,
                                            germline = gp, species = "sp",
                                            holdout_v = hold, max_rows = 16)])
    tool2 = CallableAnnotator("echo2") do seqs, ids, germline
        CallRecord[CallRecord(ids[i], seqs[i], "", "", "") for i in eachindex(ids)]
    end
    suite = suite_from_manifest(man, [tool2]; timing = TimingSpec(warmup = 0, repeats = 0))
    @test length(suite.panels) >= 3  # sim full/train/held + airr splits
    result = run_suite(suite; mode = DiagnosticMode(max_sequences = 4), store = nothing)
    @test result.mode == "diagnostic"
    @test !isempty(result.metrics)
end

@testset "SwiftIGAnnotator missing binary" begin
    a = SwiftIGAnnotator(; bin = "/nonexistent/swiftig_bin_xyz")
    @test_throws ErrorException annotate(a, ["ACGT"], ["1"], GermlinePaths(; v = VFA, d = DFA, j = JFA))
end

@testset "IgBLASTAnnotator without IgBLAST" begin
    @test_throws ErrorException IgBLASTAnnotator()
end

@testset "PanelCache freeze" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    src = SimSource(; id = "cache_toy", germline = gp, species = "sp", n = 6, seed = 9)
    panel = SimGoldPanel(src, :full; n = 6)
    cdir = mktempdir()
    cache = PanelCache(cdir)
    d1 = load_panel_cached(panel, cache)
    d2 = load_panel_cached(panel, cache)
    @test d1.sequences == d2.sequences
    @test d1.ids == d2.ids
    @test !isnothing(d1.gold) && d1.gold == d2.gold
end
