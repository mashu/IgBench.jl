# Compare IgBLAST vs swig-cli RIAT-MP and AER on an IgSim gold panel,
# plus a D-ablation copy of the same panel (90% of gold D remnants spliced out).
#
# Annotated grid only: IgBLAST + aux, swig-cli v0.37.2 with prepare-reference.
# Assignment-only and v0.34.0 are omitted (same calls, slower).
#
# IgSim is GitHub `main` (https://github.com/mashu/IgSim.jl). This script
# only overrides `indel_rate` (default 0.3% of reads).
#
# Tools (each is annotated once per panel and timed):
#   igblast-aux      IgBLAST + J/CDR3 aux
#   swig-riat-prep   v0.37.2 --assigner riat_mp --prepared-reference
#   swig-aer-prep    v0.37.2 --assigner aer --prepared-reference
#
#   export IGBENCH_V=... IGBENCH_D=... IGBENCH_J=... IGBENCH_AUX=...
#   export IGBENCH_N=100000 IGBENCH_THREADS=8
#   export IGBENCH_DROP_D=0.9
#   export IGBENCH_INDEL_P=0.003
#   julia --project=. examples/compare_igblast_swig.jl

using IgBench
using IgSim
using Downloads
using JSON

const DEFAULT_V = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_V.fasta"
const DEFAULT_D = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_D.fasta"
const DEFAULT_J = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_J.fasta"
const DEFAULT_AUX = "/home/mateusz/Mac-IgBLAST/data/rhesus_monkey_gl.aux"

const ROOT = dirname(@__DIR__)
const SWIG_VERSION = get(ENV, "SWIG_NEW", "0.37.2")
const N = parse(Int, get(ENV, "IGBENCH_N", "100000"))
const THREADS = parse(Int, get(ENV, "IGBENCH_THREADS", "8"))
const DROP_D = parse(Float64, get(ENV, "IGBENCH_DROP_D", "0.9"))
const INDEL_P = parse(Float64, get(ENV, "IGBENCH_INDEL_P", "0.003"))
const ORGANISM = get(ENV, "IGBENCH_ORGANISM", "rhesus_monkey")
const SPECIES = get(ENV, "IGBENCH_SPECIES", "rhesus")
const AUX = get(ENV, "IGBENCH_AUX", DEFAULT_AUX)
const OUT = get(ENV, "IGBENCH_OUT", joinpath(ROOT, "runs", "full", "igblast_swig"))

function comparison_params()
    (0 <= INDEL_P <= 1) || error("IGBENCH_INDEL_P must be in [0, 1], got $INDEL_P")
    default_indel = IgSim.train_params().indel_rate
    IgSim.train_params(; indel_rate = IgSim.GatedRate(INDEL_P, default_indel.rate))
end

function swig_release_asset(version::AbstractString)
    ver = lstrip(String(version), 'v')
    if Sys.islinux() && Sys.ARCH === :x86_64
        cpu = isfile("/proc/cpuinfo") ? read("/proc/cpuinfo", String) : ""
        return occursin("avx2", cpu) ? "swig-cli-linux-x64-modern" : "swig-cli-linux-x64"
    elseif Sys.islinux() && Sys.ARCH === :aarch64
        return "swig-cli-linux-arm64"
    elseif Sys.isapple() && Sys.ARCH === :aarch64
        return "swig-cli-macos-arm64"
    elseif Sys.isapple()
        return "swig-cli-macos-x64"
    elseif Sys.iswindows()
        return "swig-cli-windows-x64.exe"
    end
    error("no swig-cli build for $(Sys.KERNEL) $(Sys.ARCH)")
end

function ensure_swig_cli(version::AbstractString)
    ver = lstrip(String(version), 'v')
    asset = swig_release_asset(ver)
    dest = joinpath(ROOT, "deps", "swig", "v$ver", asset)
    if isfile(dest)
        println("using $dest")
        return dest
    end
    mkpath(dirname(dest))
    url = "https://github.com/MurrellGroup/swig/releases/download/v$ver/$asset"
    println("downloading $url")
    Downloads.download(url, dest)
    chmod(dest, 0o755)
    dest
end

function ensure_prepared_reference(bin, gp::GermlinePaths, organism::AbstractString)
    prefix = joinpath(ROOT, "deps", "swig", "refs", organism)
    manifest = prefix * ".swig-reference.json"
    if isfile(manifest)
        println("using prepared reference $manifest")
        return manifest
    end
    mkpath(dirname(prefix))
    cmd = `$bin prepare-reference -germline_db_V $(gp.v) -germline_db_D $(gp.d) -germline_db_J $(gp.j) -organism $organism -ig_seqtype Ig --out-prefix $prefix`
    println("prepare-reference → $prefix")
    success(cmd) || error("prepare-reference failed")
    isfile(manifest) || error("prepare-reference did not write $manifest")
    manifest
end

function merge_timing_rows(saved, new_rows)
    by = Dict{Tuple{String,String},Any}()
    for t in saved
        by[(String(t["panel"]), String(t["tool"]))] = t
    end
    for row in new_rows
        by[(String(row["panel"]), String(row["tool"]))] = row
    end
    keys_sorted = sort!(collect(keys(by)))
    Dict{String,Any}[by[k] for k in keys_sorted]
end

function restore_saved_timing!(store, result, saved)
    out = merge_timing_rows(saved, result.timing)
    IgBench.write_timing!(store, out)
    gallery_path = joinpath(store.root, "span_gallery.json")
    gallery = isfile(gallery_path) ? JSON.parsefile(gallery_path) : Any[]
    manifest = JSON.parsefile(joinpath(store.root, "manifest.json"))
    payload = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "suite" => result.name,
        "mode" => result.mode,
        "step" => result.step,
        "tags" => result.tags,
        "tools" => manifest["tools"],
        "panels" => JSON.parsefile(joinpath(store.root, "panels.json")),
        "metrics" => result.metrics,
        "timing" => out,
        "span_gallery" => gallery,
    )
    IgBench.write_report_if_full(FullReportMode(), store, payload)
    out
end

function print_allele_by_panel(metrics)
    rows = [r for r in metrics if r["metric"] == "allele" && String(r["ref"]) == ":gold"]
    isempty(rows) && return
    panels = unique(String(r["panel"]) for r in rows)
    println()
    println("allele vs gold")
    println(rpad("tool", 18), join((rpad(p, 36) for p in panels), " "), "  ΔV (pp)")
    tools = unique(String(r["pred"]) for r in rows)
    for tool in tools
        print(rpad(tool, 18))
        vs = Float64[]
        for p in panels
            found = false
            for r in rows
                String(r["panel"]) == p && String(r["pred"]) == tool || continue
                cell = "V=$(round(100 * Float64(r["v"]); digits=2)) D=$(round(100 * Float64(r["d"]); digits=2)) J=$(round(100 * Float64(r["j"]); digits=2))"
                print(rpad(cell, 36))
                push!(vs, Float64(r["v"]))
                found = true
                break
            end
            found || print(rpad("—", 36))
        end
        if length(vs) == 2
            print("  ", round(100 * (vs[2] - vs[1]); digits = 2))
        end
        println()
    end
end

const BIN = ensure_swig_cli(SWIG_VERSION)
println("swig-cli v$SWIG_VERSION → $BIN")
println("assigners riat_mp, aer (prepared-reference)")
println("V ", get(ENV, "IGBENCH_V", DEFAULT_V))
println("D ", get(ENV, "IGBENCH_D", DEFAULT_D))
println("J ", get(ENV, "IGBENCH_J", DEFAULT_J))
println("aux ", AUX)
println("threads ", THREADS)
println("drop_d_frac ", DROP_D)
println("indel_p ", INDEL_P, " (", round(100 * INDEL_P; digits = 2), "% of reads)")
println("IgSim ", pkgdir(IgSim))

gp = GermlinePaths(;
    v = get(ENV, "IGBENCH_V", DEFAULT_V),
    d = get(ENV, "IGBENCH_D", DEFAULT_D),
    j = get(ENV, "IGBENCH_J", DEFAULT_J),
)

man = DatasetManifest("igblast_swig";
    sim = [SimSource(; id = "sim_igsim_indelp$(INDEL_P)", db_label = "assign", germline = gp,
                     species = SPECIES, n = N, seed = 1,
                     params_factory = comparison_params)],
    airr = AirrSource[],
)

function swig_prep(name, assigner, prep)
    SwiftIGAnnotator(;
        name,
        bin = BIN,
        threads = THREADS,
        extra_args = String["--vdj", "--assigner", assigner,
                            "--prepared-reference", prep],
    )
end

PREP = ensure_prepared_reference(BIN, gp, ORGANISM)
println("prepared-reference ", PREP)
isfile(AUX) || error("IGBENCH_AUX is not a file: $AUX")

current_tools = Set(["igblast-aux", "swig-riat-prep", "swig-aer-prep"])
saved_timing = isfile(joinpath(OUT, "timing.json")) ?
    filter(t -> String(t["tool"]) in current_tools,
           JSON.parsefile(joinpath(OUT, "timing.json"))) : Any[]

tools = AbstractAnnotator[
    IgBLASTAnnotator(; name = "igblast-aux", organism_param = ORGANISM,
                     aux = AUX, num_threads = THREADS),
    swig_prep("swig-riat-prep", "riat_mp", PREP),
    swig_prep("swig-aer-prep", "aer", PREP),
]

store = DirectoryRunStore(OUT)
suite = suite_from_manifest(man, tools; drop_d_frac = DROP_D)
result = run_suite(suite;
                   mode = FullReportMode(),
                   store,
                   cache = PanelCache(joinpath(OUT, "panel_cache_igsim_indelp$(INDEL_P)")),
                   reuse_predictions = true)
timing = restore_saved_timing!(store, result, saved_timing)
println("wrote $(length(result.metrics)) metric rows → $(result.store_path)")
println("tools: ", join(tool_name.(tools), ", "))
print_allele_by_panel(result.metrics)
for t in timing
    println(t["panel"], "  ", t["tool"], "  ", round(Float64(t["wall_s"]); digits = 3), " s  ",
            round(Float64(t["seq_per_s"]); digits = 1), " seq/s")
end
println(joinpath(result.store_path, "report.html"))
