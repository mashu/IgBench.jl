# Compare IgBLAST vs swig-cli RIAT-MP, AER, and AER+igblast_balanced
# on IgSim gold panels (IgM and IgG SHM priors).
#
# Annotated grid: IgBLAST + aux, swig-cli v0.37.2 with prepare-reference.
# All tools are timed and scored on the same reads (allele + spans).
#
# IgSim is GitHub `main` (https://github.com/mashu/IgSim.jl). This script
# overrides `indel_rate` (default 0.3% of reads) and SHM isotype
# (`shm_igm` / `shm_igg`). No constant-region sequence is simulated.
#
# Tools (each is annotated once and timed):
#   igblast-aux                 IgBLAST + J/CDR3 aux
#   swig-riat-prep              v0.37.2 --assigner riat_mp --prepared-reference
#   swig-aer-prep               v0.37.2 --assigner aer --prepared-reference
#   swig-aer-igblast-balanced   same AER path plus --calling-profile igblast_balanced
#
#   export IGBENCH_V=... IGBENCH_D=... IGBENCH_J=... IGBENCH_AUX=...
#   export IGBENCH_N=20000 IGBENCH_THREADS=8
#   export IGBENCH_INDEL_P=0.003
#   export IGBENCH_ISOTYPES=igm,igg
#   julia --project=. examples/compare_igblast_swig.jl

using IgBench
using IgSim
using Downloads

const DEFAULT_V = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_V.fasta"
const DEFAULT_D = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_D.fasta"
const DEFAULT_J = "/home/mateusz/Mac-IgBLAST/data/Macaca_mulatta_J.fasta"
const DEFAULT_AUX = "/home/mateusz/Mac-IgBLAST/data/rhesus_monkey_gl.aux"

const ROOT = dirname(@__DIR__)
const SWIG_VERSION = get(ENV, "SWIG_NEW", "0.37.2")
const N = parse(Int, get(ENV, "IGBENCH_N", "20000"))
const THREADS = parse(Int, get(ENV, "IGBENCH_THREADS", "8"))
const INDEL_P = parse(Float64, get(ENV, "IGBENCH_INDEL_P", "0.003"))
const ORGANISM = get(ENV, "IGBENCH_ORGANISM", "rhesus_monkey")
const SPECIES = get(ENV, "IGBENCH_SPECIES", "rhesus")
const AUX = get(ENV, "IGBENCH_AUX", DEFAULT_AUX)
const OUT = get(ENV, "IGBENCH_OUT", joinpath(ROOT, "runs", "full", "igblast_swig"))

function parse_isotypes(raw::AbstractString)
    out = String[]
    for part in split(raw, ',')
        s = lowercase(strip(String(part)))
        isempty(s) && continue
        (s == "igm" || s == "igg") ||
            error("IGBENCH_ISOTYPES entries must be igm or igg, got $(part)")
        s in out || push!(out, s)
    end
    isempty(out) && error("IGBENCH_ISOTYPES is empty")
    out
end

const ISOTYPES = parse_isotypes(get(ENV, "IGBENCH_ISOTYPES", "igm,igg"))

function panel_source_id(isotype::AbstractString)
    isotype == "igm" ?
        "sim_igsim_indelp$(INDEL_P)_n$(N)" :
        "sim_igsim_$(isotype)_indelp$(INDEL_P)_n$(N)"
end

function comparison_params(isotype::AbstractString)
    (0 <= INDEL_P <= 1) || error("IGBENCH_INDEL_P must be in [0, 1], got $INDEL_P")
    default_indel = IgSim.train_params().indel_rate
    shm = isotype == "igg" ? IgSim.shm_igg() : IgSim.shm_igm()
    IgSim.train_params(; body_error_rate = shm,
                       indel_rate = IgSim.GatedRate(INDEL_P, default_indel.rate))
end

function isotype_label(isotype::AbstractString)
    isotype == "igg" ? "IgG" : "IgM"
end

function seed_panel_cache!(dest::AbstractString, src::AbstractString)
    isdir(src) || return
    mkpath(dest)
    for f in readdir(src)
        dst = joinpath(dest, f)
        isfile(dst) && continue
        srcf = joinpath(src, f)
        isfile(srcf) && cp(srcf, dst)
    end
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

function print_allele_by_panel(metrics)
    rows = [r for r in metrics if r["metric"] == "allele" && String(r["ref"]) == ":gold"]
    isempty(rows) && return
    panels = unique(String(r["panel"]) for r in rows)
    println()
    println("allele vs gold")
    println(rpad("tool", 28), join((rpad(p, 36) for p in panels), " "))
    tools = unique(String(r["pred"]) for r in rows)
    for tool in tools
        print(rpad(tool, 28))
        for p in panels
            found = false
            for r in rows
                String(r["panel"]) == p && String(r["pred"]) == tool || continue
                cell = "V=$(round(100 * Float64(r["v"]); digits=2)) D=$(round(100 * Float64(r["d"]); digits=2)) J=$(round(100 * Float64(r["j"]); digits=2))"
                print(rpad(cell, 36))
                found = true
                break
            end
            found || print(rpad("—", 36))
        end
        println()
    end
end

const BIN = ensure_swig_cli(SWIG_VERSION)
println("swig-cli v$SWIG_VERSION → $BIN")
println("assigners riat_mp, aer, aer+igblast_balanced (prepared-reference)")
println("V ", get(ENV, "IGBENCH_V", DEFAULT_V))
println("D ", get(ENV, "IGBENCH_D", DEFAULT_D))
println("J ", get(ENV, "IGBENCH_J", DEFAULT_J))
println("aux ", AUX)
println("n ", N)
println("threads ", THREADS)
println("indel_p ", INDEL_P, " (", round(100 * INDEL_P; digits = 2), "% of reads)")
println("isotypes ", join(ISOTYPES, ", "))
println("IgSim ", pkgdir(IgSim))

gp = GermlinePaths(;
    v = get(ENV, "IGBENCH_V", DEFAULT_V),
    d = get(ENV, "IGBENCH_D", DEFAULT_D),
    j = get(ENV, "IGBENCH_J", DEFAULT_J),
)

man = DatasetManifest("igblast_swig";
    sim = [let iso = iso
               SimSource(; id = panel_source_id(iso), db_label = "assign",
                         germline = gp, species = SPECIES, n = N, seed = 1,
                         label = isotype_label(iso),
                         params_factory = () -> comparison_params(iso))
           end for iso in ISOTYPES],
    airr = AirrSource[],
)

function swig_prep(name, assigner, prep, extra = String[])
    SwiftIGAnnotator(;
        name,
        bin = BIN,
        threads = THREADS,
        extra_args = String["--vdj", "--assigner", assigner, extra...,
                            "--prepared-reference", prep],
    )
end

PREP = ensure_prepared_reference(BIN, gp, ORGANISM)
println("prepared-reference ", PREP)
isfile(AUX) || error("IGBENCH_AUX is not a file: $AUX")

tools = AbstractAnnotator[
    IgBLASTAnnotator(; name = "igblast-aux", organism_param = ORGANISM,
                     aux = AUX, num_threads = THREADS),
    swig_prep("swig-riat-prep", "riat_mp", PREP),
    swig_prep("swig-aer-prep", "aer", PREP),
    swig_prep("swig-aer-igblast-balanced", "aer", PREP,
              String["--calling-profile", "igblast_balanced"]),
]

store = DirectoryRunStore(OUT)
cache_dir = joinpath(OUT, "panel_cache")
seed_panel_cache!(cache_dir, joinpath(OUT, "panel_cache_sim_igsim_indelp$(INDEL_P)_n$(N)"))
suite = suite_from_manifest(man, tools)
result = run_suite(suite;
                   mode = FullReportMode(),
                   store,
                   cache = PanelCache(cache_dir),
                   reuse_predictions = true)
println("wrote $(length(result.metrics)) metric rows → $(result.store_path)")
println("tools: ", join(tool_name.(tools), ", "))
print_allele_by_panel(result.metrics)
for t in result.timing
    println(t["panel"], "  ", t["tool"], "  ", round(Float64(t["wall_s"]); digits = 3), " s  ",
            round(Float64(t["seq_per_s"]); digits = 1), " seq/s")
end
println(joinpath(result.store_path, "report.html"))
