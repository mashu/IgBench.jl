# panel.jl — Loaded evaluation cohorts (sim gold or real AIRR).

"""Loaded panel ready for annotate + metrics."""
struct PanelData
    id::String
    species::String
    germline::GermlinePaths
    sequences::Vector{String}
    ids::Vector{String}
    gold::Union{Nothing,Vector{CallRecord}}
    file_gold::Union{Nothing,Vector{CallRecord}}
    meta::Dict{String,Any}
end

"""Something that can be loaded into [`PanelData`](@ref)."""
abstract type AbstractPanel end

panel_id(p::AbstractPanel) = p.id

"""Convert IgSim labeled read → [`CallRecord`](@ref)."""
function callrecord_from_labeled(r::IgSim.LabeledRead, id::AbstractString)
    sp = r.spans
    CallRecord(
        String(id),
        r.sequence,
        r.v_call,
        ismissing(r.d_call) ? "" : String(r.d_call),
        r.j_call,
        Span(sp.v.start, sp.v.stop),
        Span(sp.d.start, sp.d.stop),
        Span(sp.j.start, sp.j.stop),
    )
end

"""
Simulated panel with causal gold labels.

`gallery` ∈ `:full` (any V from full DB), `:train` (held removed), `:held`
(only held V via `IgSim.HoldoutVGenerator`).
"""
struct SimGoldPanel <: AbstractPanel
    id::String
    source::SimSource
    gallery::Symbol
    n::Int
end

function SimGoldPanel(source::SimSource, gallery::Symbol; n::Union{Nothing,Integer} = nothing)
    gallery in (:full, :train, :held) ||
        error("gallery must be :full, :train, or :held; got $gallery")
    id = "$(source.id)_$(gallery)"
    SimGoldPanel(id, source, gallery, isnothing(n) ? source.n : Int(n))
end

function load_panel(p::SimGoldPanel)
    src = p.source
    gp = src.germline
    full_db = IgSim.load_germline(; v = gp.v, d = gp.d, j = gp.j)
    params = src.params_factory()
    train_db, held_db = IgSim.holdout_alleles(full_db, src.holdout_v, String[], String[])
    if p.gallery === :full
        gen = IgSim.ReadGenerator(full_db; params)
    elseif p.gallery === :train
        gen = IgSim.ReadGenerator(train_db; params)
    else
        isempty(src.holdout_v) && error("SimGoldPanel gallery=:held needs non-empty holdout_v")
        gen = IgSim.HoldoutVGenerator(full_db, src.holdout_v; params)
    end
    rng = MersenneTwister(src.seed + (p.gallery === :full ? 0 : p.gallery === :train ? 1 : 2))
    reads = Vector{IgSim.LabeledRead}(undef, p.n)
    @inbounds for i in 1:p.n
        reads[i] = gen(rng)
    end
    ids = ["$(p.id)_$i" for i in 1:length(reads)]
    gold = CallRecord[callrecord_from_labeled(reads[i], ids[i]) for i in eachindex(reads)]
    seqs = String[r.sequence for r in reads]
    PanelData(p.id, src.species, gp, seqs, ids, gold, nothing,
              Dict{String,Any}(
                  "kind" => "sim",
                  "gallery" => String(p.gallery),
                  "species" => src.species,
                  "n" => length(seqs),
                  "seed" => src.seed,
                  "source_id" => src.id,
              ))
end

"""
Real AIRR panel.

`membership` ∈ `:all`, `:seen` (V not in holdout), `:held` (V in holdout).
Requires `holdout_v` on the source for `:seen`/`:held`.
"""
struct AirrPanel <: AbstractPanel
    id::String
    source::AirrSource
    membership::Symbol
    n::Union{Nothing,Int}
end

function AirrPanel(source::AirrSource, membership::Symbol = :all;
                   n::Union{Nothing,Integer} = nothing)
    membership in (:all, :seen, :held) ||
        error("membership must be :all, :seen, or :held; got $membership")
    id = membership === :all ? source.id : "$(source.id)_$(membership)"
    AirrPanel(id, source, membership, isnothing(n) ? source.max_rows : Int(n))
end

function v_in_holdout(v_call::AbstractString, holdout::AbstractSet{String})
    for tok in parse_allele_calls(v_call)
        normalize_allele(tok) in holdout && return true
        tok in holdout && return true
    end
    false
end

function load_panel(p::AirrPanel)
    src = p.source
    rows = read_airr_calls(src.path; max_rows = nothing, skip_nonproductive = false)
    hold = Set(normalize_allele(x) for x in src.holdout_v)
    union!(hold, Set(src.holdout_v))
    if p.membership === :seen
        isempty(src.holdout_v) && error("AirrPanel membership=:seen needs holdout_v")
        rows = CallRecord[r for r in rows if !v_in_holdout(r.v_call, hold)]
    elseif p.membership === :held
        isempty(src.holdout_v) && error("AirrPanel membership=:held needs holdout_v")
        rows = CallRecord[r for r in rows if v_in_holdout(r.v_call, hold)]
    end
    if !isnothing(p.n) && length(rows) > p.n
        rows = rows[1:p.n]
    end
    ids = String[r.sequence_id for r in rows]
    seqs = String[r.sequence for r in rows]
    file_gold = rows  # calls from file
    # gold for metrics against file labels
    PanelData(p.id, src.species, src.germline, seqs, ids, nothing, file_gold,
              Dict{String,Any}(
                  "kind" => "airr",
                  "membership" => String(p.membership),
                  "species" => src.species,
                  "organism_param" => src.organism_param,
                  "n" => length(seqs),
                  "path" => src.path,
                  "source_id" => src.id,
              ))
end

"""Cap requested panel size for the run mode."""
function panel_for_mode(p::SimGoldPanel, mode::RunMode)
    n_eff = effective_n(mode, p.n)
    n_eff == p.n ? p : SimGoldPanel(p.source, p.gallery; n = n_eff)
end

function panel_for_mode(p::AirrPanel, mode::DiagnosticMode)
    n_eff = isnothing(p.n) ? mode.max_sequences : effective_n(mode, p.n)
    AirrPanel(p.source, p.membership; n = n_eff)
end

panel_for_mode(p::AirrPanel, ::FullReportMode) = p
panel_for_mode(p::AbstractPanel, ::RunMode) = p

panel_germline(p::SimGoldPanel) = p.source.germline
panel_germline(p::AirrPanel) = p.source.germline
panel_species(p::SimGoldPanel) = p.source.species
panel_species(p::AirrPanel) = p.source.species

"""
Disk cache of frozen panel sequences + labels (AIRR TSV.gz + meta JSON).

Ensures mid-training diagnostics see the same N reads every step. Pass to
[`run_suite`](@ref) as `cache=PanelCache("cache/demo")`.
"""
struct PanelCache
    root::String
end

PanelCache(root::AbstractString) = PanelCache(String(root))

function cache_paths(cache::PanelCache, panel_id::AbstractString)
    safe = replace(String(panel_id), r"[^\w\.-]" => "_")
    (airr = joinpath(cache.root, safe * ".airr.tsv.gz"),
     meta = joinpath(cache.root, safe * ".meta.json"))
end

function panel_data_from_rows(p::AbstractPanel, rows::Vector{CallRecord},
                              meta::Dict{String,Any})
    ids = String[r.sequence_id for r in rows]
    seqs = String[r.sequence for r in rows]
    kind = String(get(meta, "kind", ""))
    if kind == "sim"
        PanelData(panel_id(p), panel_species(p), panel_germline(p),
                  seqs, ids, rows, nothing, meta)
    else
        PanelData(panel_id(p), panel_species(p), panel_germline(p),
                  seqs, ids, nothing, rows, meta)
    end
end

function rows_for_cache(data::PanelData)
    !isnothing(data.gold) && return data.gold
    !isnothing(data.file_gold) && return data.file_gold
    CallRecord[CallRecord(data.ids[i], data.sequences[i], "", "", "")
               for i in eachindex(data.ids)]
end

"""Load panel, reading/writing [`PanelCache`](@ref) when provided."""
load_panel_cached(p::AbstractPanel, ::Nothing) = load_panel(p)

function load_panel_cached(p::AbstractPanel, cache::PanelCache)
    paths = cache_paths(cache, panel_id(p))
    if isfile(paths.airr) && isfile(paths.meta)
        rows = read_airr_calls(paths.airr)
        meta = Dict{String,Any}(String(k) => v for (k, v) in JSON.parsefile(paths.meta))
        return panel_data_from_rows(p, rows, meta)
    end
    data = load_panel(p)
    mkpath(cache.root)
    write_airr_calls(paths.airr, rows_for_cache(data))
    open(paths.meta, "w") do io
        JSON.print(io, data.meta, 2)
    end
    data
end
