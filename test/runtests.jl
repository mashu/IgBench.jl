using Test
using IgBench
import IgSim
using Random

const FIX = joinpath(@__DIR__, "fixtures")
const VFA = joinpath(FIX, "V.fasta")
const DFA = joinpath(FIX, "D.fasta")
const JFA = joinpath(FIX, "J.fasta")

@testset "allele score (comma only)" begin
    @test parse_allele_calls("IGHV1-69*01/IGHV1-69D*01") ==
          ["IGHV1-69*01/IGHV1-69D*01"]
    @test parse_allele_calls("A,B") == ["A", "B"]
    @test parse_allele_calls("A, B, C") == ["A", "B", "C"]
    @test primary_allele_call("IGHV1-46*01,IGHV1-46*03") == "IGHV1-46*01"
    @test primary_allele_call("IGHV1-69*01/IGHV1-69D*01") ==
          "IGHV1-69*01/IGHV1-69D*01"
    @test allele_score("A", "A") == 1
    @test allele_score("A,B", "A") ≈ 1 / 2
    @test allele_score("A,B,C", "A") ≈ 1 / 3
    @test allele_score("C", "A") == 0
    @test allele_score("A/B", "A") == 0
    @test allele_score("A/B", "A/B") == 1
    @test allele_score("A,B", "A,B") == 1
    @test allele_score("", "IGHV1-46*01") == 0
    @test allele_score("IGHV1-69*01/IGHV1-69D*01",
                       "IGHV1-69*01/IGHV1-69D*01") == 1
    @test call_field_empty("") && call_field_empty("NA") && call_field_empty(".")
    names = metric_name.(default_metrics())
    @test names == ["allele", "span_iou", "span_exact", "span_start", "span_stop"]
end

@testset "spans" begin
    @test span_iou(Span(1, 10), Span(5, 15)) ≈ 6 / 15
    @test span_iou(EMPTY_SPAN, Span(1, 5)) == 0.0
    @test isnan(span_iou(Span(1, 5), EMPTY_SPAN))
    gold = Span(27, 314)
    pred = Span(31, 314)
    @test span_iou(pred, gold) ≈ 284 / 288
    @test span_exact(pred, gold) == 0
    @test span_start(pred, gold) == 0
    @test span_stop(pred, gold) == 1
    @test span_exact(EMPTY_SPAN, gold) == 0
    @test isnan(span_exact(pred, EMPTY_SPAN))
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

@testset "AlleleAccuracy + span metrics" begin
    gold = [
        CallRecord("a", "AAA", "V1*01", "D1*01", "J1*01", Span(1, 1), Span(2, 2), Span(3, 3)),
        CallRecord("b", "BBB", "V2*01", "", "J2*01", Span(1, 1), EMPTY_SPAN, Span(2, 2)),
    ]
    pred = [
        CallRecord("a", "AAA", "V1*01", "D1*01", "J1*01", Span(1, 1), Span(2, 2), Span(3, 3)),
        CallRecord("b", "BBB", "V2*01", "", "J9*01", Span(1, 1), EMPTY_SPAN, Span(2, 2)),
    ]
    m = evaluate(AlleleAccuracy(), pred, gold)
    @test m.v == 1.0
    @test m.j == 0.5
    @test m.d == 1.0
    @test m.n == 2
    @test m.d_n == 1
    @test m.v_n == 2 && m.j_n == 2

    gold2 = [
        CallRecord("a", "AAA", "A", "D", "J"),
        CallRecord("b", "BBB", "", "X", "J"),
    ]
    pred2 = [
        CallRecord("a", "AAA", "A,B", "D", "J"),
        CallRecord("b", "BBB", "Z", "", "J"),
    ]
    f = evaluate(AlleleAccuracy(), pred2, gold2)
    @test f.v ≈ 0.5 && f.v_n == 1
    @test f.d == 0.5 && f.d_n == 2
    @test f.j == 1.0 && f.j_n == 2
    cm = call_metrics(["A,B"], [""], ["J"], ["A"], ["X"], ["J"])
    @test cm.v ≈ 0.5 && cm.d == 0.0 && cm.j == 1.0
    @test cm isa MetricValue

    gold_span = [
        CallRecord("a", "AAA", "V", "D", "J", Span(27, 314), Span(10, 20), Span(30, 40)),
        CallRecord("b", "BBB", "V", "D", "J", Span(1, 10), Span(11, 15), Span(16, 20)),
    ]
    pred_span = [
        CallRecord("a", "AAA", "V", "D", "J", Span(31, 314), Span(10, 20), Span(30, 40)),
        CallRecord("b", "BBB", "V", "", "J", EMPTY_SPAN, EMPTY_SPAN, Span(16, 20)),
    ]
    iou = evaluate(SpanIoU(), pred_span, gold_span)
    @test iou.v ≈ (284 / 288 + 0.0) / 2
    @test iou.d == 0.5
    @test iou.j == 1.0
    ex = evaluate(SpanExact(), pred_span, gold_span)
    @test ex.v == 0.0
    @test ex.d == 0.5
    @test ex.j == 1.0
    st = evaluate(SpanStart(), pred_span, gold_span)
    @test st.v == 0.0
    @test st.d == 0.5
    sp = evaluate(SpanStop(), pred_span, gold_span)
    @test sp.v == 0.5
    recs = call_records(["A"], ["D"], ["J"]; prefix = "x")
    @test length(recs) == 1 && recs[1].sequence_id == "x1" && recs[1].v_call == "A"
end

@testset "SimGoldPanel + run_suite Fake" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    db = IgSim.load_germline(; v = VFA, d = DFA, j = JFA)
    hold = String[db.v[1].name]
    src = SimSource(; id = "toy", db_label = "fix", germline = gp, species = "toy_species",
                    holdout_v = hold, n = 8, seed = 1)
    panel = SimGoldPanel(src, :all; n = 8)
    @test occursin("assign=fix", panel.id)
    @test occursin("sim=all", panel.id)
    @test !occursin("held", panel.id) && !occursin("seen", panel.id) &&
          !occursin("closed", panel.id)
    data = load_panel(panel)
    @test length(data.sequences) == 8
    @test !isnothing(data.gold)
    @test data.meta["sim_set"] == "all"
    @test data.meta["assign_v"] == VFA

    closed = SimGoldPanel(src, :minus_holdout; n = 8, assign = gp, assign_label = "train")
    @test occursin("assign=train", closed.id)
    @test occursin("sim=minus_holdout", closed.id)
    cdata = load_panel(closed)
    @test cdata.germline.v == VFA
    @test cdata.meta["sim_set"] == "minus_holdout"
    @test cdata.meta["db_label"] == "train"
    @test cdata.meta["sim_v"] == VFA

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
    @test any(r -> r["metric"] == "allele" && r["v"] == 1.0, result.metrics)
    @test any(r -> r["metric"] == "span_exact", result.metrics)
    @test isfile(joinpath(outdir, "metrics.jsonl"))
    @test isfile(joinpath(outdir, "metrics.json"))
    @test isfile(joinpath(outdir, "summary.md"))
    @test isfile(joinpath(outdir, "manifest.json"))
    @test !isempty(result.timing)
    @test haskey(result.timing[1], "wall_s")
    @test haskey(result.timing[1], "n_sequences")
    @test haskey(result.timing[1], "seq_per_s")
end

@testset "suite_from_manifest + diagnostic" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    db = IgSim.load_germline(; v = VFA, d = DFA, j = JFA)
    hold = String[db.v[1].name]
    src = SimSource(; id = "toy", db_label = "fix", germline = gp, species = "sp",
                    holdout_v = hold, n = 16, seed = 2)
    data = load_panel(SimGoldPanel(src, :all; n = 16))
    airr_path = joinpath(tempdir(), "igbench_mini.airr.tsv")
    write_airr_calls(airr_path, data.gold)
    man = DatasetManifest("mini";
                          sim = [src],
                          airr = [AirrSource(; id = "real_toy", db_label = "fix",
                                            path = airr_path,
                                            germline = gp, species = "sp",
                                            max_rows = 16)])
    tool2 = CallableAnnotator("echo2") do seqs, ids, germline
        CallRecord[CallRecord(ids[i], seqs[i], "", "", "") for i in eachindex(ids)]
    end
    suite = suite_from_manifest(man, [tool2]; timing = TimingSpec(warmup = 0, repeats = 0))
    @test length(suite.panels) == 4  # sim all/minus_holdout/holdout_only + one AIRR
    @test all(p -> occursin("assign=", p.id), suite.panels)
    result = run_suite(suite; mode = DiagnosticMode(max_sequences = 4), store = nothing)
    @test result.mode == "diagnostic"
    @test !isempty(result.metrics)
end

@testset "AIRR gold validator" begin
    gp = GermlinePaths(; v = VFA, d = DFA, j = JFA)
    db = IgSim.load_germline(; v = VFA, d = DFA, j = JFA)
    ok = [
        CallRecord("s1", "ACGT", db.v[1].name, db.d[1].name, db.j[1].name),
    ]
    names = germline_allele_names(gp)
    validate_gold_calls!(ok, names, "ok")
    bad = [
        CallRecord("s1", "ACGT", "NOT-IN-DB*01", db.d[1].name, db.j[1].name),
    ]
    @test_throws ErrorException validate_gold_calls!(bad, names, "bad")
    airr_path = joinpath(tempdir(), "igbench_stale.airr.tsv")
    write_airr_calls(airr_path, bad)
    src = AirrSource(; id = "stale", db_label = "fix", path = airr_path,
                     germline = gp, species = "sp")
    @test_throws ErrorException load_panel(AirrPanel(src))
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
    src = SimSource(; id = "cache_toy", db_label = "fix", germline = gp, species = "sp", n = 6, seed = 9)
    panel = SimGoldPanel(src, :all; n = 6)
    cdir = mktempdir()
    cache = PanelCache(cdir)
    d1 = load_panel_cached(panel, cache)
    d2 = load_panel_cached(panel, cache)
    @test d1.sequences == d2.sequences
    @test d1.ids == d2.ids
    @test !isnothing(d1.gold) && d1.gold == d2.gold
end
