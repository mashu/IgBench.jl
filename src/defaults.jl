# defaults.jl — Expand DatasetManifest into panels + default CompareSpecs.

"""
Build default compares: each tool vs `:gold`, and tool↔tool pairs.
"""
function default_compares(tools;
                          vs_gold::Bool = true,
                          pairs::Bool = true)
    names = [tool_name(t) for t in tools]
    comps = CompareSpec[]
    metrics = default_metrics()
    agree = agreement_metrics()
    for n in names
        vs_gold && push!(comps, CompareSpec(n, ":gold", metrics))
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

Sim sources → `sim_set ∈ (:all, :minus_holdout, :holdout_only)`
(`minus_holdout` / `holdout_only` only if holdout is non-empty).
AIRR sources → one panel each (assign FASTA on the source).
"""
function suite_from_manifest(manifest::DatasetManifest,
                             tools;
                             name::AbstractString = manifest.name,
                             timing::TimingSpec = TimingSpec(),
                             sim_sets = SIM_SETS,
                             drop_d_frac::Real = 0)
    panels = AbstractPanel[]
    for src in manifest.sim
        for s in sim_sets
            if s in (:minus_holdout, :holdout_only) && isempty(src.holdout_v)
                continue
            end
            push!(panels, SimGoldPanel(src, s))
        end
    end
    for src in manifest.airr
        push!(panels, AirrPanel(src))
    end
    frac = Float64(drop_d_frac)
    if frac > 0
        extras = AbstractPanel[]
        for p in panels
            append_drop_d!(extras, p, frac)
        end
        append!(panels, extras)
    end
    comps = default_compares(tools)
    BenchSuite(name; panels, tools = collect(tools), compares = comps, timing)
end
