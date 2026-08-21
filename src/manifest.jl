# manifest.jl — Species-agnostic dataset description (paths only; nothing bundled).

"""
Real AIRR library: path + assign germlines + optional gold germlines.

`germline` is the FASTA given to tools (assign_db). `gold_germline` defaults to
the same paths; gold call tokens must be names in that FASTA.
"""
struct AirrSource
    id::String
    db_label::String
    path::String
    germline::GermlinePaths
    gold_germline::GermlinePaths
    species::String
    organism_param::String
    max_rows::Union{Nothing,Int}
end

function AirrSource(;
                    id::AbstractString,
                    path::AbstractString,
                    germline::GermlinePaths,
                    species::AbstractString,
                    db_label::AbstractString = id,
                    gold_germline = nothing,
                    organism_param::AbstractString = species,
                    max_rows = nothing)
    gold = isnothing(gold_germline) ? germline : gold_germline
    AirrSource(String(id), String(db_label), String(path), germline, gold,
               String(species), String(organism_param),
               isnothing(max_rows) ? nothing : Int(max_rows))
end

"""
Simulated source: assign FASTA + holdout list + sample size.

`germline` is the FASTA given to tools (assign_db) and the parent DB IgSim
samples from. `sim_set` on the panel chooses which alleles are drawn.
"""
struct SimSource{F}
    id::String
    db_label::String
    germline::GermlinePaths
    species::String
    holdout_v::Vector{String}
    n::Int
    seed::Int
    label::String
    params_factory::F
end

function SimSource(;
                   id::AbstractString,
                   germline::GermlinePaths,
                   species::AbstractString,
                   db_label::AbstractString = id,
                   holdout_v = String[],
                   n::Integer = 512,
                   seed::Integer = 1,
                   label::AbstractString = "",
                   params_factory = IgSim.train_params)
    SimSource(String(id), String(db_label), germline, String(species),
              String[String(x) for x in holdout_v],
              Int(n), Int(seed), String(label), params_factory)
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

function germline_to_dict(g::GermlinePaths)
    Dict{String,Any}("v" => g.v, "d" => g.d, "j" => g.j)
end

function manifest_to_dict(m::DatasetManifest)
    Dict{String,Any}(
        "name" => m.name,
        "sim" => [Dict{String,Any}(
            "id" => s.id,
            "db_label" => s.db_label,
            "species" => s.species,
            "n" => s.n,
            "seed" => s.seed,
            "label" => s.label,
            "holdout_v" => s.holdout_v,
            "germline" => germline_to_dict(s.germline),
        ) for s in m.sim],
        "airr" => [Dict{String,Any}(
            "id" => a.id,
            "db_label" => a.db_label,
            "path" => a.path,
            "species" => a.species,
            "organism_param" => a.organism_param,
            "max_rows" => a.max_rows,
            "germline" => germline_to_dict(a.germline),
            "gold_germline" => germline_to_dict(a.gold_germline),
        ) for a in m.airr],
    )
end
