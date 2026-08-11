# defaults.jl — Expand DatasetManifest into panels + default CompareSpecs.

"""
Build default compares: each tool vs `:gold` / `:file_gold`, and tool↔tool pairs.
"""
function default_compares(tools;
                          vs_gold::Bool = true,
                          vs_file_gold::Bool = true,
                          pairs::Bool = true)
    names = [tool_name(t) for t in tools]
    comps = CompareSpec[]
    metrics = default_metrics()
    agree = agreement_metrics()
    for n in names
        vs_gold && push!(comps, CompareSpec(n, ":gold", metrics))
        vs_file_gold && push!(comps, CompareSpec(n, ":file_gold", metrics))
    end
    if pairs
        for i in 1:length(names)
            for j in i+1:length(names)
                push!(comps, CompareSpec(names[i], names[j], agree))
            end
        end
    end
    comps
end

"""
Expand a [`DatasetManifest`](@ref) into a [`BenchSuite`](@ref).

Sim sources → `gallery ∈ (:full,:train,:held)` (held/train only if holdout non-empty).
AIRR sources → `membership ∈ (:all,:seen,:held)` similarly.
"""
function suite_from_manifest(manifest::DatasetManifest,
                             tools;
                             name::AbstractString = manifest.name,
                             timing::TimingSpec = TimingSpec(),
                             galleries = (:full, :train, :held),
                             memberships = (:all, :seen, :held))
    panels = AbstractPanel[]
    for src in manifest.sim
        for g in galleries
            if g in (:train, :held) && isempty(src.holdout_v)
                continue
            end
            push!(panels, SimGoldPanel(src, g))
        end
    end
    for src in manifest.airr
        for m in memberships
            if m in (:seen, :held) && isempty(src.holdout_v)
                continue
            end
            push!(panels, AirrPanel(src, m))
        end
    end
    comps = default_compares(tools)
    BenchSuite(name; panels, tools = collect(tools), compares = comps, timing)
end
