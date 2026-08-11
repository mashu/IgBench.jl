# manifest.jl — Species-agnostic dataset description (paths only; nothing bundled).

"""Real AIRR cohort: path + germlines + free species tag."""
struct AirrSource
    id::String
    path::String
    germline::GermlinePaths
    species::String
    organism_param::String
    holdout_v::Vector{String}
    max_rows::Union{Nothing,Int}
end

function AirrSource(;
                    id::AbstractString,
                    path::AbstractString,
                    germline::GermlinePaths,
                    species::AbstractString,
                    organism_param::AbstractString = species,
                    holdout_v = String[],
                    max_rows = nothing)
    AirrSource(String(id), String(path), germline, String(species),
               String(organism_param), String[String(x) for x in holdout_v],
               isnothing(max_rows) ? nothing : Int(max_rows))
end

"""Simulated cohort: IgSim germlines + holdout list + sample size."""
struct SimSource{F}
    id::String
    germline::GermlinePaths
    species::String
    holdout_v::Vector{String}
    n::Int
    seed::Int
    params_factory::F
end

function SimSource(;
                   id::AbstractString,
                   germline::GermlinePaths,
                   species::AbstractString,
                   holdout_v = String[],
                   n::Integer = 512,
                   seed::Integer = 1,
                   params_factory = IgSim.train_params)
    SimSource(String(id), germline, String(species),
              String[String(x) for x in holdout_v],
              Int(n), Int(seed), params_factory)
end

"""Named collection of sim + AIRR sources."""
struct DatasetManifest
    name::String
    sim::Vector{<:SimSource}
    airr::Vector{AirrSource}
end

DatasetManifest(name::AbstractString;
                sim = SimSource[],
                airr = AirrSource[]) =
    DatasetManifest(String(name), collect(sim), collect(airr))

function manifest_to_dict(m::DatasetManifest)
    Dict{String,Any}(
        "name" => m.name,
        "sim" => [Dict{String,Any}(
            "id" => s.id,
            "species" => s.species,
            "n" => s.n,
            "seed" => s.seed,
            "holdout_v" => s.holdout_v,
            "germline" => Dict("v" => s.germline.v, "d" => s.germline.d, "j" => s.germline.j),
        ) for s in m.sim],
        "airr" => [Dict{String,Any}(
            "id" => a.id,
            "path" => a.path,
            "species" => a.species,
            "organism_param" => a.organism_param,
            "holdout_v" => a.holdout_v,
            "max_rows" => a.max_rows,
            "germline" => Dict("v" => a.germline.v, "d" => a.germline.d, "j" => a.germline.j),
        ) for a in m.airr],
    )
end
