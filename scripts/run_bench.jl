#!/usr/bin/env julia
# Standalone entry: parse args → load manifest → run_suite.
#
# Usage:
#   julia --project=. scripts/run_bench.jl --manifest path/to/manifest.jl \
#         --mode full|diagnostic --out runs/full/demo [--tools igblast,swiftig]

using IgBench
using Dates

# Package extension activates only after IgBLAST is loaded.
IgBench.load_igblast_extension!()

function parse_args(args)
    manifest = ""
    mode = "full"
    out = ""
    tools = ["igblast"]
    organism = "human"
    i = 1
    while i <= length(args)
        if args[i] == "--manifest" && i < length(args)
            manifest = args[i + 1]; i += 2
        elseif args[i] == "--mode" && i < length(args)
            mode = args[i + 1]; i += 2
        elseif args[i] == "--out" && i < length(args)
            out = args[i + 1]; i += 2
        elseif args[i] == "--tools" && i < length(args)
            tools = split(args[i + 1], ','); i += 2
        elseif args[i] == "--organism" && i < length(args)
            organism = args[i + 1]; i += 2
        else
            error("unknown arg: $(args[i])")
        end
    end
    isempty(manifest) && error("--manifest is required")
    isempty(out) && (out = joinpath("runs", mode, string(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))))
    (; manifest, mode, out, tools, organism)
end

function build_tools(names, organism)
    out = AbstractAnnotator[]
    for n in names
        n = String(strip(n))
        if n == "igblast"
            push!(out, IgBLASTAnnotator(; organism_param = organism))
        elseif n == "swiftig"
            push!(out, SwiftIGAnnotator())
        else
            error("unknown tool '$n' (built-ins: igblast, swiftig; use library API for CallableAnnotator)")
        end
    end
    out
end

function main(args)
    opt = parse_args(args)
    # Manifest file should return a DatasetManifest, e.g.:
    #   using IgBench
    #   DatasetManifest(...)
    man = include(abspath(opt.manifest))
    man isa DatasetManifest || error("manifest must evaluate to DatasetManifest")
    tools = build_tools(opt.tools, opt.organism)
    suite = suite_from_manifest(man, tools)
    mode = opt.mode == "diagnostic" ? DiagnosticMode() : FullReportMode()
    store = DirectoryRunStore(opt.out)
    result = run_suite(suite; mode, store)
    println("wrote $(length(result.metrics)) metric rows → $(result.store_path)")
    opt.mode == "full" && println(joinpath(result.store_path, "report.html"))
    result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
