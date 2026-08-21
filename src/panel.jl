# panel.jl — Loaded evaluation sets (sim gold or real AIRR).

const SIM_SETS = (:all, :minus_holdout, :holdout_only)

"""Loaded panel ready for annotate + metrics."""
struct PanelData
    id::String
    species::String
    germline::GermlinePaths
    sequences::Vector{String}
    ids::Vector{String}
    gold::Union{Nothing,Vector{CallRecord}}
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

function sim_panel_id(source_id::AbstractString, assign_label::AbstractString, sim_set::Symbol)
    "$(source_id)__assign=$(assign_label)__sim=$(sim_set)"
end

function sim_panel_id(source::SimSource, sim_set::Symbol)
    sim_panel_id(source.id, source.db_label, sim_set)
end

function airr_panel_id(source::AirrSource)
    "$(source.id)__assign=$(source.db_label)"
end

"""Allele names in a germline FASTA (IgSim header parsing)."""
function germline_allele_names(gp::GermlinePaths)
    db = IgSim.load_germline(; v = gp.v, d = gp.d, j = gp.j)
    names = Set{String}()
    for a in db.v
        push!(names, a.name)
    end
    for a in db.d
        push!(names, a.name)
    end
    for a in db.j
        push!(names, a.name)
    end
    names
end

"""
Every non-empty comma-token of V/D/J calls must be a FASTA name.

Hard error if the TSV was annotated with a different database.
"""
function validate_gold_calls!(rows::AbstractVector{CallRecord}, names::AbstractSet{<:AbstractString},
                              panel_id::AbstractString)
    for r in rows
        for (field, call) in (("v_call", r.v_call), ("d_call", r.d_call), ("j_call", r.j_call))
            call_field_empty(call) && continue
            for tok in parse_allele_calls(call)
                tok in names && continue
                error("panel $(panel_id): $field token '$tok' is not in the gold germline FASTA; " *
                      "this TSV was not annotated with this database")
            end
        end
    end
    rows
end

"""
Simulated panel with causal gold labels.

`sim_set` ∈ `:all` (any allele from **source** FASTA), `:minus_holdout`
(holdout V removed), `:holdout_only` (only holdout V).

`assign` is the FASTA given to tools (default: `source.germline`). Pass a
smaller FASTA to score “sim from A, assign to B”.
"""
struct SimGoldPanel <: AbstractPanel
    id::String
    source::SimSource
    sim_set::Symbol
    n::Int
    assign::GermlinePaths
    assign_label::String
end

function SimGoldPanel(source::SimSource, sim_set::Symbol;
                      n::Union{Nothing,Integer} = nothing,
                      assign = nothing,
                      assign_label::Union{Nothing,AbstractString} = nothing)
    sim_set in SIM_SETS ||
        error("sim_set must be one of $SIM_SETS; got $sim_set")
    gp = isnothing(assign) ? source.germline : assign
    label = isnothing(assign_label) ? source.db_label : String(assign_label)
    SimGoldPanel(sim_panel_id(source.id, label, sim_set), source, sim_set,
                 isnothing(n) ? source.n : Int(n), gp, label)
end

function sim_set_seed_offset(sim_set::Symbol)
    sim_set === :all && return 0
    sim_set === :minus_holdout && return 1
    2
end

function load_panel(p::SimGoldPanel)
    src = p.source
    gp = src.germline
    full_db = IgSim.load_germline(; v = gp.v, d = gp.d, j = gp.j)
    params = src.params_factory()
    train_db, _held_db = IgSim.holdout_alleles(full_db, src.holdout_v, String[], String[])
    if p.sim_set === :all
        gen = IgSim.ReadGenerator(full_db; params)
    elseif p.sim_set === :minus_holdout
        isempty(src.holdout_v) && error("SimGoldPanel sim_set=:minus_holdout needs non-empty holdout_v")
        gen = IgSim.ReadGenerator(train_db; params)
    else
        isempty(src.holdout_v) && error("SimGoldPanel sim_set=:holdout_only needs non-empty holdout_v")
        gen = IgSim.HoldoutVGenerator(full_db, src.holdout_v; params)
    end
    rng = MersenneTwister(src.seed + sim_set_seed_offset(p.sim_set))
    reads = Vector{IgSim.LabeledRead}(undef, p.n)
    @inbounds for i in 1:p.n
        reads[i] = gen(rng)
    end
    ids = ["$(p.id)_$i" for i in 1:length(reads)]
    gold = CallRecord[callrecord_from_labeled(reads[i], ids[i]) for i in eachindex(reads)]
    seqs = String[r.sequence for r in reads]
    agp = p.assign
    meta = Dict{String,Any}(
        "kind" => "sim",
        "sim_set" => String(p.sim_set),
        "db_label" => p.assign_label,
        "assign_v" => agp.v,
        "assign_d" => agp.d,
        "assign_j" => agp.j,
        "sim_v" => gp.v,
        "sim_d" => gp.d,
        "sim_j" => gp.j,
        "species" => src.species,
        "n" => length(seqs),
        "seed" => src.seed,
        "source_id" => src.id,
        "holdout_v" => src.holdout_v,
    )
    !isempty(src.label) && (meta["label"] = src.label)
    PanelData(p.id, src.species, agp, seqs, ids, gold, meta)
end

"""
Real AIRR panel. Gold is the file's V/D/J calls after validating that every
token is a name in `gold_germline` (default: assign FASTA).
"""
struct AirrPanel <: AbstractPanel
    id::String
    source::AirrSource
    n::Union{Nothing,Int}
end

function AirrPanel(source::AirrSource; n::Union{Nothing,Integer} = nothing)
    AirrPanel(airr_panel_id(source), source, isnothing(n) ? source.max_rows : Int(n))
end

function load_panel(p::AirrPanel)
    src = p.source
    rows = read_airr_calls(src.path; max_rows = nothing, skip_nonproductive = false)
    names = germline_allele_names(src.gold_germline)
    validate_gold_calls!(rows, names, p.id)
    if !isnothing(p.n) && length(rows) > p.n
        rows = rows[1:p.n]
    end
    ids = String[r.sequence_id for r in rows]
    seqs = String[r.sequence for r in rows]
    gp = src.germline
    gold_gp = src.gold_germline
    PanelData(p.id, src.species, gp, seqs, ids, rows,
              Dict{String,Any}(
                  "kind" => "airr",
                  "db_label" => src.db_label,
                  "assign_v" => gp.v,
                  "assign_d" => gp.d,
                  "assign_j" => gp.j,
                  "gold_v" => gold_gp.v,
                  "gold_d" => gold_gp.d,
                  "gold_j" => gold_gp.j,
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
    n_eff == p.n ? p : SimGoldPanel(p.source, p.sim_set; n = n_eff,
                                    assign = p.assign, assign_label = p.assign_label)
end

function panel_for_mode(p::AirrPanel, mode::DiagnosticMode)
    n_eff = isnothing(p.n) ? mode.max_sequences : effective_n(mode, p.n)
    AirrPanel(p.source; n = n_eff)
end

panel_for_mode(p::AirrPanel, ::FullReportMode) = p
panel_for_mode(p::AbstractPanel, ::RunMode) = p

panel_germline(p::SimGoldPanel) = p.assign
panel_germline(p::AirrPanel) = p.source.germline
panel_species(p::SimGoldPanel) = p.source.species
panel_species(p::AirrPanel) = p.source.species

"""Short report name; empty means use the panel id."""
panel_report_label(::AbstractPanel) = ""
panel_report_label(p::SimGoldPanel) = p.source.label

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
    PanelData(panel_id(p), panel_species(p), panel_germline(p),
              seqs, ids, rows, meta)
end

function rows_for_cache(data::PanelData)
    !isnothing(data.gold) && return data.gold
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
