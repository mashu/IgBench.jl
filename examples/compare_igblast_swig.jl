# Compare IgBLAST vs swig-cli on an IgSim gold panel, with annotation modes
# timed as separate tools (aux / CDR3 annotation changes wall-clock).
#
# Downloads swig-cli v0.37.2 into deps/swig/ on first run (override with
# SWIFTIG_BIN or SWIG_VERSION). IgBLAST.jl must already be in the environment.
#
# Tools (each is annotated once and timed):
#   igblast-off      IgBLAST, no -auxiliary_data
#   igblast-aux      IgBLAST + NCBI/custom aux (J/CDR3)
#   swig-off         swig-cli --vdj (assignment only; no FWR/CDR/junction)
#   swig-aux         same aux file via -auxiliary_data (IgBLAST-style J metadata)
#   swig-annots      --prepared-reference (Swig metadata from prepare-reference)
#
# swig-cli prepare-reference is run once (cached under deps/swig/refs/).
# --prepared-reference cannot mix with germline FASTAs or -auxiliary_data.
#
#   export IGBENCH_AUX=/path/to/other.aux          # optional override
#   export IGBENCH_V=... IGBENCH_D=... IGBENCH_J=...
#   export IGBENCH_N=100000 IGBENCH_THREADS=8
#   julia --project=. examples/compare_igblast_swig.jl

using IgBench
using Downloads
using JSON

# Machine defaults (this host). Override with IGBENCH_* / SWIFTIG_BIN / SWIG_VERSION.
const DEFAULT_V = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_V.fasta"
const DEFAULT_D = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_D.fasta"
const DEFAULT_J = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_J.fasta"
const DEFAULT_AUX = "/home/mateusz/Mac-IgBLAST/data/rhesus_monkey_gl.aux"

const ROOT = dirname(@__DIR__)
const SWIG_VERSION = get(ENV, "SWIG_VERSION", "0.37.2")
const N = parse(Int, get(ENV, "IGBENCH_N", "100000"))
const THREADS = parse(Int, get(ENV, "IGBENCH_THREADS", "8"))
const ORGANISM = get(ENV, "IGBENCH_ORGANISM", "rhesus_monkey")
const SPECIES = get(ENV, "IGBENCH_SPECIES", "rhesus")
const AUX = get(ENV, "IGBENCH_AUX", DEFAULT_AUX)
const ASSIGNER = get(ENV, "IGBENCH_ASSIGNER", "riat_mp")
const OUT = get(ENV, "IGBENCH_OUT", joinpath(ROOT, "runs", "full", "igblast_swig"))

function swig_release_asset()
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

"""Download pinned swig-cli into deps/swig/ unless SWIFTIG_BIN is set."""
function ensure_swig_cli(version::AbstractString)
    override = get(ENV, "SWIFTIG_BIN", "")
    !isempty(override) && return override
    asset = swig_release_asset()
    dest = joinpath(ROOT, "deps", "swig", "v$version", asset)
    if isfile(dest)
        println("using $dest")
        return dest
    end
    mkpath(dirname(dest))
    url = "https://github.com/MurrellGroup/swig/releases/download/v$version/$asset"
    println("downloading $url")
    Downloads.download(url, dest)
    chmod(dest, 0o755)
    dest
end

"""Infer Swig V/D/J metadata once; reuse via --prepared-reference."""
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

function prediction_path(root, tool)
    dir = joinpath(root, "predictions")
    isdir(dir) || return ""
    suffix = "__$(replace(String(tool), r"[^\w\.-]" => "_")).airr.tsv.gz"
    for f in readdir(dir)
        endswith(f, suffix) && return joinpath(dir, f)
    end
    ""
end

function replay_annotator(name, path)
    println("reusing IgBLAST calls $path")
    rows = read_airr_calls(path; max_rows = nothing, skip_nonproductive = false)
    by_id = Dict(r.sequence_id => r for r in rows)
    CallableAnnotator(name) do seqs, ids, germline
        n = length(ids)
        out = Vector{CallRecord}(undef, n)
        for i in 1:n
            r = get(by_id, String(ids[i]), nothing)
            out[i] = isnothing(r) ? CallRecord(ids[i], seqs[i], "", "", "") : r
        end
        out
    end
end

function igblast_tool(name; aux = "")
    path = prediction_path(OUT, name)
    !isempty(path) && return replay_annotator(name, path)
    if isempty(aux)
        return IgBLASTAnnotator(; name, organism_param = ORGANISM, num_threads = THREADS)
    end
    IgBLASTAnnotator(; name, organism_param = ORGANISM, aux, num_threads = THREADS)
end

function restore_igblast_timing!(store, result, saved)
    isempty(saved) && return result.timing
    by = Dict{String,Any}(String(t["tool"]) => t for t in saved)
    out = Dict{String,Any}[]
    for row in result.timing
        name = String(row["tool"])
        if startswith(name, "igblast") && haskey(by, name)
            merged = copy(row)
            old = by[name]
            merged["wall_s"] = old["wall_s"]
            merged["seq_per_s"] = old["seq_per_s"]
            merged["n_sequences"] = old["n_sequences"]
            push!(out, merged)
        else
            push!(out, row)
        end
    end
    IgBench.write_timing!(store, out)
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
    )
    IgBench.write_report_if_full(FullReportMode(), store, payload)
    out
end

const BIN = ensure_swig_cli(SWIG_VERSION)
println("swig-cli v$SWIG_VERSION → $BIN")
println("V ", get(ENV, "IGBENCH_V", DEFAULT_V))
println("D ", get(ENV, "IGBENCH_D", DEFAULT_D))
println("J ", get(ENV, "IGBENCH_J", DEFAULT_J))
println("aux ", AUX)
println("threads ", THREADS)

gp = GermlinePaths(;
    v = get(ENV, "IGBENCH_V", DEFAULT_V),
    d = get(ENV, "IGBENCH_D", DEFAULT_D),
    j = get(ENV, "IGBENCH_J", DEFAULT_J),
)

PREP = ensure_prepared_reference(BIN, gp, ORGANISM)
println("prepared-reference ", PREP)

man = DatasetManifest("igblast_swig";
    sim = [SimSource(; id = "sim", db_label = "assign", germline = gp,
                     species = SPECIES, n = N, seed = 1)],
    airr = AirrSource[],
)

swig(name, extra) = SwiftIGAnnotator(;
    name,
    bin = BIN,
    threads = THREADS,
    extra_args = String["--vdj", "--assigner", ASSIGNER,
                        "-organism", ORGANISM, "-ig_seqtype", "Ig", extra...],
)

saved_timing = isfile(joinpath(OUT, "timing.json")) ?
    JSON.parsefile(joinpath(OUT, "timing.json")) : Any[]

tools = AbstractAnnotator[
    igblast_tool("igblast-off"),
    swig("swig-off", String[]),
    swig("swig-annots", ["--prepared-reference", PREP]),
]

if !isempty(AUX)
    isfile(AUX) || error("IGBENCH_AUX is not a file: $AUX")
    push!(tools, igblast_tool("igblast-aux"; aux = AUX))
    push!(tools, swig("swig-aux", ["-auxiliary_data", AUX]))
else
    println("IGBENCH_AUX unset — skipping igblast-aux and swig-aux")
end

store = DirectoryRunStore(OUT)
suite = suite_from_manifest(man, tools)
result = run_suite(suite;
                   mode = FullReportMode(),
                   store,
                   cache = PanelCache(joinpath(OUT, "panel_cache")))
timing = restore_igblast_timing!(store, result, saved_timing)
println("wrote $(length(result.metrics)) metric rows → $(result.store_path)")
println("tools: ", join(tool_name.(tools), ", "))
for t in timing
    println(t["tool"], "  ", round(Float64(t["wall_s"]); digits = 3), " s  ",
            round(Float64(t["seq_per_s"]); digits = 1), " seq/s")
end
println(joinpath(result.store_path, "report.html"))
