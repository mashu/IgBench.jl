# drop_d.jl — Ablation panel: splice gold D remnants out of existing reads.

"""Shift a span after deleting the closed interval `del` from the sequence."""
function span_after_delete(s::Span, del::Span)
    isempty(s) && return EMPTY_SPAN
    n = length(del)
    s.stop < del.start && return s
    s.start > del.stop && return Span(s.start - n, s.stop - n)
    left = s.start < del.start
    right = s.stop > del.stop
    left && !right && return Span(s.start, del.start - 1)
    right && !left && return Span(del.stop + 1 - n, s.stop - n)
    left && right && return Span(s.start, s.stop - n)
    EMPTY_SPAN
end

"""Remove the gold D interval from one record; keep V/J alleles, clear D gold."""
function drop_d_from_record(r::CallRecord)
    del = r.d_span
    isempty(del) && return CallRecord(r.sequence_id, r.sequence, r.v_call, "", r.j_call,
                                      r.v_span, EMPTY_SPAN, r.j_span)
    seq = r.sequence
    (del.start < 1 || del.stop > length(seq)) &&
        error("d_span $(del.start)-$(del.stop) outside sequence length $(length(seq)) for $(r.sequence_id)")
    newseq = seq[1:del.start - 1] * seq[del.stop + 1:end]
    CallRecord(r.sequence_id, newseq, r.v_call, "", r.j_call,
               span_after_delete(r.v_span, del), EMPTY_SPAN,
               span_after_delete(r.j_span, del))
end

"""
Copy a loaded sim panel, splicing D out of `drop_frac` of reads that have a
gold D span. V/J gold alleles stay; D gold is cleared on dropped reads.
"""
function drop_d_from_panel(data::PanelData, drop_frac::Float64, seed::Int;
                          id::AbstractString = data.id)
    isnothing(data.gold) && error("drop-D panel $(data.id) has no gold")
    (0 <= drop_frac <= 1) || error("drop_frac must be in [0, 1], got $drop_frac")
    gold = data.gold
    rng = MersenneTwister(seed)
    n_with = 0
    n_drop = 0
    new_gold = Vector{CallRecord}(undef, length(gold))
    for i in eachindex(gold)
        r = gold[i]
        if isempty(r.d_span)
            new_gold[i] = r
            continue
        end
        n_with += 1
        if rand(rng) < drop_frac
            n_drop += 1
            new_gold[i] = drop_d_from_record(r)
        else
            new_gold[i] = r
        end
    end
    seqs = String[r.sequence for r in new_gold]
    meta = Dict{String,Any}(String(k) => v for (k, v) in data.meta)
    meta["kind"] = "sim_drop_d"
    meta["drop_d_frac"] = drop_frac
    meta["drop_d_seed"] = seed
    meta["d_present_before"] = n_with
    meta["d_dropped"] = n_drop
    meta["parent_id"] = data.id
    meta["n"] = length(seqs)
    PanelData(String(id), data.species, data.germline, seqs, data.ids, new_gold, meta)
end

"""Derived panel: same reads as `parent`, with a fraction of gold D remnants excised."""
struct DropDPanel{P<:AbstractPanel} <: AbstractPanel
    id::String
    parent::P
    drop_frac::Float64
    seed::Int
end

function DropDPanel(parent::AbstractPanel; drop_frac::Real = 0.9, seed::Integer = 1)
    frac = Float64(drop_frac)
    (0 <= frac <= 1) || error("drop_frac must be in [0, 1], got $frac")
    DropDPanel{typeof(parent)}(panel_id(parent) * "__drop_d=$(frac)",
                               parent, frac, Int(seed))
end

panel_id(p::DropDPanel) = p.id

load_panel(p::DropDPanel) =
    drop_d_from_panel(load_panel(p.parent), p.drop_frac, p.seed; id = p.id)

function load_panel_cached(p::DropDPanel, cache::PanelCache)
    paths = cache_paths(cache, panel_id(p))
    if isfile(paths.airr) && isfile(paths.meta)
        rows = read_airr_calls(paths.airr)
        meta = Dict{String,Any}(String(k) => v for (k, v) in JSON.parsefile(paths.meta))
        return panel_data_from_rows(p, rows, meta)
    end
    parent_data = load_panel_cached(p.parent, cache)
    data = drop_d_from_panel(parent_data, p.drop_frac, p.seed; id = p.id)
    mkpath(cache.root)
    write_airr_calls(paths.airr, rows_for_cache(data))
    open(paths.meta, "w") do io
        JSON.print(io, data.meta, 2)
    end
    data
end

function panel_for_mode(p::DropDPanel, mode::RunMode)
    DropDPanel(panel_for_mode(p.parent, mode); drop_frac = p.drop_frac, seed = p.seed)
end

panel_germline(p::DropDPanel) = panel_germline(p.parent)
panel_species(p::DropDPanel) = panel_species(p.parent)
panel_report_label(p::DropDPanel) = panel_report_label(p.parent)

function append_drop_d!(out::Vector{AbstractPanel}, p::SimGoldPanel, frac::Float64)
    push!(out, DropDPanel(p; drop_frac = frac, seed = p.source.seed))
    nothing
end

append_drop_d!(::Vector{AbstractPanel}, ::AbstractPanel, ::Float64) = nothing
