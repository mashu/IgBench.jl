# span_gallery.jl — Sample start/stop disagreements onto a shared MSA.

const SPAN_SAMPLE_N = 10

record_span(r::CallRecord, ::Val{:v}) = r.v_span
record_span(r::CallRecord, ::Val{:d}) = r.d_span
record_span(r::CallRecord, ::Val{:j}) = r.j_span

record_call(r::CallRecord, ::Val{:v}) = r.v_call
record_call(r::CallRecord, ::Val{:d}) = r.d_call
record_call(r::CallRecord, ::Val{:j}) = r.j_call

kind_score(::Val{:start}, pred::Span, gold::Span) = span_start(pred, gold)
kind_score(::Val{:stop}, pred::Span, gold::Span) = span_stop(pred, gold)
kind_score(::Val{:iou}, pred::Span, gold::Span) = span_iou(pred, gold)

locus_key(::Val{:v}) = "v"
locus_key(::Val{:d}) = "d"
locus_key(::Val{:j}) = "j"
kind_key(::Val{:start}) = "start"
kind_key(::Val{:stop}) = "stop"
kind_key(::Val{:iou}) = "iou"

function kind_mismatch(kind, pred::Span, gold::Span)
    s = kind_score(kind, pred, gold)
    !isnan(s) && s < 1
end

"""Allele name → DNA string from assign FASTA (V, D, and J merged)."""
function germline_sequence_index(gp::GermlinePaths)
    db = IgSim.load_germline(; v = gp.v, d = gp.d, j = gp.j)
    idx = Dict{String,String}()
    for a in db.v
        idx[a.name] = a.sequence
    end
    for a in db.d
        idx[a.name] = a.sequence
    end
    for a in db.j
        idx[a.name] = a.sequence
    end
    idx
end

function lookup_germline(idx::AbstractDict{String,String}, call::AbstractString)
    name = primary_allele_call(call)
    isempty(name) && return "", ""
    haskey(idx, name) || return name, ""
    name, idx[name]
end

"""
Semi-global alignment: consume the whole query; free gaps at subject ends.

Returns `(query_aln, subject_aln, match_line)`.
"""
function semiglobal_align(query::AbstractString, subject::AbstractString)
    q = uppercase(String(query))
    s = uppercase(String(subject))
    n = ncodeunits(q)
    m = ncodeunits(s)
    (n == 0 || m == 0) && return q, s, ""

    F = Matrix{Int}(undef, n + 1, m + 1)
    ptr = Matrix{UInt8}(undef, n + 1, m + 1)
    F[1, 1] = 0
    ptr[1, 1] = 0x00
    for j in 1:m
        F[1, j + 1] = 0
        ptr[1, j + 1] = 0x02
    end
    for i in 1:n
        F[i + 1, 1] = i * ALIGN_GAP
        ptr[i + 1, 1] = 0x01
    end
    for i in 1:n
        qi = codeunit(q, i)
        for j in 1:m
            sj = codeunit(s, j)
            sub = qi == sj ? ALIGN_MATCH : ALIGN_MISMATCH
            diag = F[i, j] + sub
            up = F[i, j + 1] + ALIGN_GAP
            left = F[i + 1, j] + ALIGN_GAP
            if diag >= up && diag >= left
                F[i + 1, j + 1] = diag
                ptr[i + 1, j + 1] = 0x00
            elseif up >= left
                F[i + 1, j + 1] = up
                ptr[i + 1, j + 1] = 0x01
            else
                F[i + 1, j + 1] = left
                ptr[i + 1, j + 1] = 0x02
            end
        end
    end
    best_j = 0
    best = typemin(Int)
    for j in 0:m
        sc = F[n + 1, j + 1]
        if sc >= best
            best = sc
            best_j = j
        end
    end
    traceback_align(q, s, ptr, n, best_j)
end

function traceback_align(q::AbstractString, s::AbstractString, ptr::Matrix{UInt8},
                         i::Int, j_end::Int)
    m = ncodeunits(s)
    j = j_end
    qa = Char[]
    sa = Char[]
    while i > 0
        p = ptr[i + 1, j + 1]
        if p == 0x00
            push!(qa, Char(codeunit(q, i)))
            push!(sa, Char(codeunit(s, j)))
            i -= 1
            j -= 1
        elseif p == 0x01
            push!(qa, Char(codeunit(q, i)))
            push!(sa, '-')
            i -= 1
        else
            push!(qa, '-')
            push!(sa, Char(codeunit(s, j)))
            j -= 1
        end
    end
    reverse!(qa)
    reverse!(sa)
    prefix = j == 0 ? "" : SubString(s, 1, j)
    suffix = j_end >= m ? "" : SubString(s, j_end + 1, m)
    qstr = ("-" ^ ncodeunits(prefix)) * String(qa) * ("-" ^ ncodeunits(suffix))
    sstr = String(prefix) * String(sa) * String(suffix)
    qstr, sstr, match_line(qstr, sstr)
end

function sample_disagree_indices(idxs::AbstractVector{Int}, n::Integer, rng::AbstractRNG)
    k = min(Int(n), length(idxs))
    k == 0 && return Int[]
    length(idxs) == k && return collect(idxs)
    collect(idxs[randperm(rng, length(idxs))[1:k]])
end

"""IgBLAST vs SWIG families; other names are their own family."""
function tool_family(name::AbstractString)
    n = lowercase(String(name))
    startswith(n, "igblast") && return "igblast"
    startswith(n, "swig") && return "swig"
    n
end

span_coord(::Val{:start}, sp::Span) = isempty(sp) ? nothing : sp.start
span_coord(::Val{:stop}, sp::Span) = isempty(sp) ? nothing : sp.stop
span_coord(::Val{:iou}, sp::Span) = isempty(sp) ? nothing : (sp.start, sp.stop)

function span_delta(::Val{:start}, pred::Span, gold::Span)
    (isempty(pred) || isempty(gold)) && return nothing
    pred.start - gold.start
end

function span_delta(::Val{:stop}, pred::Span, gold::Span)
    (isempty(pred) || isempty(gold)) && return nothing
    pred.stop - gold.stop
end

function span_nt(seq::AbstractString, sp::Span)
    isempty(sp) && return ""
    n = ncodeunits(seq)
    (sp.start < 1 || sp.stop > n) && return ""
    String(SubString(seq, sp.start, sp.stop))
end

function compact_offset_hist(deltas::AbstractVector{Int}, n_missing::Integer)
    acc = Dict{Int,Int}()
    n_zero = 0
    for d in deltas
        if d == 0
            n_zero += 1
        else
            acc[d] = get(acc, d, 0) + 1
        end
    end
    ks = sort!(collect(keys(acc)))
    bins = Vector{Vector{Int}}(undef, length(ks))
    for i in eachindex(ks)
        k = ks[i]
        bins[i] = Int[k, acc[k]]
    end
    Dict{String,Any}(
        "n" => length(deltas) + Int(n_missing),
        "zeros" => n_zero,
        "missing" => Int(n_missing),
        "bins" => bins,
    )
end

function other_family_shares(preds, tool_order, focus, i, locus, kind, coord)
    fam = tool_family(focus)
    for t in tool_order
        String(t) == focus && continue
        tool_family(t) == fam && continue
        haskey(preds, t) || continue
        span_coord(kind, record_span(preds[t][i], locus)) == coord && return true
    end
    false
end

function sample_unique_shared(unique_ix, shared_ix, n_sample, rng)
    n_want = Int(n_sample)
    n_u = min(length(unique_ix), cld(n_want, 2))
    n_s = min(length(shared_ix), n_want - n_u)
    n_u = min(length(unique_ix), n_want - n_s)
    vcat(sample_disagree_indices(unique_ix, n_u, rng),
         sample_disagree_indices(shared_ix, n_s, rng))
end

function pair_frac(mat, n_scored)
    nt = size(mat, 1)
    Any[Any[a == b ? nothing : (n_scored == 0 ? nothing : mat[a, b] / n_scored) for b in 1:nt]
        for a in 1:nt]
end
const CONTRAST_SAMPLE_TOP = 6
const CONTRAST_SAMPLE_MIN_FRAC = 0.01
const CONTRAST_BIN_TOP = 3

"""Absolute edge error vs gold. Missing spans are worse than any finite offset."""
function edge_abs_error(kind, pred::Span, gold::Span)
    d = span_delta(kind, pred, gold)
    d === nothing && return typemax(Int)
    abs(Int(d))
end

function contrast_sample_ok(n_win, n_scored, rank)
    n_win > 0 || return false
    rank <= CONTRAST_SAMPLE_TOP && return true
    n_scored > 0 && n_win >= CONTRAST_SAMPLE_MIN_FRAC * n_scored
end

offset_bin_key(::Val{:start}, pred::Span, gold::Span) = delta_bin_key(span_delta(Val(:start), pred, gold))
offset_bin_key(::Val{:stop}, pred::Span, gold::Span) = delta_bin_key(span_delta(Val(:stop), pred, gold))
function offset_bin_key(::Val{:iou}, pred::Span, gold::Span)
    delta_bin_key(span_delta(Val(:start), pred, gold)) * "," *
        delta_bin_key(span_delta(Val(:stop), pred, gold))
end

delta_bin_key(::Nothing) = "miss"
delta_bin_key(d::Integer) = string(d)

function format_bin_label(::Val{:start}, key::AbstractString)
    key == "miss" ? "missing start" : "Δstart $(signed_delta(key))"
end
function format_bin_label(::Val{:stop}, key::AbstractString)
    key == "miss" ? "missing stop" : "Δstop $(signed_delta(key))"
end
function format_bin_label(::Val{:iou}, key::AbstractString)
    parts = split(String(key), ',')
    length(parts) == 2 || return String(key)
    "Δstart $(signed_delta(parts[1])), Δstop $(signed_delta(parts[2]))"
end

function signed_delta(s::AbstractString)
    s == "miss" && return "miss"
    n = tryparse(Int, s)
    n === nothing && return s
    n > 0 ? "+$n" : string(n)
end

function contrast_share(err_winner)
    err_winner == 0 ? "winner_exact" : "winner_closer"
end

"""MSA examples for reads where `winner` is strictly closer to gold than `loser`."""
function build_contrast_pack(preds, gold, idx, names, ixs, winner, loser, n_exact, n_closer_wrong,
                             locus, kind, n_sample, rng; store_samples::Bool = true)
    winner = String(winner)
    loser = String(loser)
    bin_packs = build_bin_packs(preds, gold, idx, names, ixs, winner, loser, locus, kind,
                                n_sample, rng; store_samples)
    peak_bin, peak_n = loser_peak_bin(bin_packs, loser)
    peak_samples = Any[]
    for p in bin_packs
        p["tool"] == loser && p["bin"] == peak_bin && (peak_samples = p["samples"])
    end
    Dict{String,Any}(
        "winner" => winner,
        "loser" => loser,
        "n" => length(ixs),
        "n_exact" => n_exact,
        "n_closer_wrong" => n_closer_wrong,
        "n_sample" => length(peak_samples),
        "peak" => peak_bin == "" ? "" : format_bin_label(kind, peak_bin),
        "peak_bin" => peak_bin,
        "peak_tool" => loser,
        "peak_n" => peak_n,
        "offsets" => contrast_offsets(preds, gold, (winner, loser), ixs, locus),
        "bin_packs" => bin_packs,
        "samples" => peak_samples,
        "loser_reads" => store_samples
            ? contrast_loser_reads(preds, gold, ixs, winner, loser, locus, kind)
            : Dict{String,Any}[],
    )
end

"""FASTA-ready reads for the column tool on this heatmap pair, tagged with Δ vs gold."""
function contrast_loser_reads(preds, gold, ixs, winner, loser, locus, kind)
    reads = Dict{String,Any}[]
    for i in ixs
        push!(reads, contrast_loser_read(gold[i], preds[loser][i], preds[winner][i], locus, kind))
    end
    reads
end

function contrast_loser_read(gold::CallRecord, pred::CallRecord, winner_pred::CallRecord,
                             locus, kind)
    seq = gold.sequence
    gspan = record_span(gold, locus)
    pspan = record_span(pred, locus)
    wspan = record_span(winner_pred, locus)
    gs, ge = span_bounds(gspan)
    ps, pe = span_bounds(pspan)
    Dict{String,Any}(
        "sequence_id" => gold.sequence_id,
        "sequence" => seq,
        "delta" => span_delta(kind, pspan, gspan),
        "winner_delta" => span_delta(kind, wspan, gspan),
        "gold_span" => gs == 0 ? Int[] : Int[gs, ge],
        "pred_span" => ps == 0 ? Int[] : Int[ps, pe],
        "length" => isempty(gspan) ? 0 : length(gspan),
        "v_call" => gold.v_call,
        "d_call" => gold.d_call,
        "j_call" => gold.j_call,
        "v_nt" => span_nt(seq, gold.v_span),
        "d_nt" => span_nt(seq, gold.d_span),
        "j_nt" => span_nt(seq, gold.j_span),
    )
end

function loser_peak_bin(bin_packs, loser)
    best_k, best_n = "", 0
    for p in bin_packs
        p["tool"] == loser || continue
        n = Int(p["n"])
        n > best_n || continue
        best_n = n
        best_k = String(p["bin"])
    end
    best_k, best_n
end

function contrast_bin_keep(nbin, pair_n, rank)
    nbin > 0 || return false
    rank <= CONTRAST_BIN_TOP && return true
    pair_n > 0 && nbin * 100 >= pair_n
end

"""Per-tool Δ bins for the heatmap pair, sampled so a hist click can open that offset."""
function build_bin_packs(preds, gold, idx, names, ixs, winner, loser, locus, kind, n_sample, rng;
                         store_samples::Bool)
    groups = Dict{Tuple{String,String},Vector{Int}}()
    for i in ixs
        gspan = record_span(gold[i], locus)
        for t in (winner, loser)
            k = offset_bin_key(kind, record_span(preds[t][i], locus), gspan)
            push!(get!(Vector{Int}, groups, (t, k)), i)
        end
    end
    pair_n = length(ixs)
    packs = Dict{String,Any}[]
    for t in (winner, loser)
        bins = Tuple{String,Vector{Int}}[(k, gix) for ((tt, k), gix) in groups if tt == t]
        sort!(bins; by = b -> (-length(b[2]), b[1]))
        for (rank, (k, gix)) in enumerate(bins)
            nbin = length(gix)
            nbin > 0 || continue
            samples = Dict{String,Any}[]
            if store_samples && contrast_bin_keep(nbin, pair_n, rank)
                for i in sample_disagree_indices(gix, n_sample, rng)
                    gspan = record_span(gold[i], locus)
                    ew = edge_abs_error(kind, record_span(preds[winner][i], locus), gspan)
                    peers = [(u, preds[u][i]) for u in names]
                    push!(samples, sample_record(gold[i], loser, preds[loser][i], peers, locus, idx;
                                                 share = contrast_share(ew)))
                end
            end
            push!(packs, Dict{String,Any}(
                "tool" => t,
                "bin" => k,
                "n" => nbin,
                "samples" => samples,
            ))
        end
    end
    packs
end

function contrast_offsets(preds, gold, names, ixs, locus)
    dstart = Dict{String,Vector{Int}}(t => Int[] for t in names)
    dstop = Dict{String,Vector{Int}}(t => Int[] for t in names)
    miss_start = Dict{String,Int}(t => 0 for t in names)
    miss_stop = Dict{String,Int}(t => 0 for t in names)
    for i in ixs
        gspan = record_span(gold[i], locus)
        for t in names
            pspan = record_span(preds[t][i], locus)
            ds = span_delta(Val(:start), pspan, gspan)
            dp = span_delta(Val(:stop), pspan, gspan)
            if isnothing(ds)
                miss_start[t] += 1
            else
                push!(dstart[t], ds)
            end
            if isnothing(dp)
                miss_stop[t] += 1
            else
                push!(dstop[t], dp)
            end
        end
    end
    Dict{String,Any}(t => Dict{String,Any}(
        "start" => compact_offset_hist(dstart[t], miss_start[t]),
        "stop" => compact_offset_hist(dstop[t], miss_stop[t]),
    ) for t in names)
end

function span_bounds(sp::Span)
    isempty(sp) ? (0, 0) : (sp.start, sp.stop)
end

function sample_record(gold::CallRecord, focus::AbstractString, pred::CallRecord,
                       peers, locus, idx::AbstractDict{String,String};
                       share::AbstractString = "")
    gspan = record_span(gold, locus)
    pspan = record_span(pred, locus)
    gcall = record_call(gold, locus)
    pcall = record_call(pred, locus)
    gname, gseq = lookup_germline(idx, gcall)
    pname, pseq = lookup_germline(idx, pcall)
    iou = span_iou(pspan, gspan)
    gs, ge = span_bounds(gspan)
    trims = Tuple{String,Int,Int,String}[("gold", gs, ge, gcall)]
    tool_spans = Dict{String,Any}()
    for (tname, rec) in peers
        ts, te = span_bounds(record_span(rec, locus))
        push!(trims, (String(tname), ts, te, record_call(rec, locus)))
        tool_spans[String(tname)] = ts == 0 ? Int[] : Int[ts, te]
    end
    Dict{String,Any}(
        "sequence_id" => gold.sequence_id,
        "gold_call" => gcall,
        "pred_call" => pcall,
        "gold_gl" => gname,
        "pred_gl" => pname,
        "gold_span" => gs == 0 ? Int[] : Int[gs, ge],
        "pred_span" => isempty(pspan) ? Int[] : Int[pspan.start, pspan.stop],
        "tool_spans" => tool_spans,
        "delta_start" => (isempty(pspan) || isempty(gspan)) ? nothing : pspan.start - gspan.start,
        "delta_stop" => (isempty(pspan) || isempty(gspan)) ? nothing : pspan.stop - gspan.stop,
        "iou" => isnan(iou) ? nothing : iou,
        "gold_gl_seq" => gseq,
        "pred_gl_seq" => isempty(pseq) ? gseq : pseq,
        "msa" => query_germline_msa(gold.sequence, gname, gseq, trims),
        "focus" => String(focus),
        "share" => share,
    )
end

"""
One gallery cell: disagreement rate + up to `n_sample` aligned examples.

`rate` is mean `(1 - score)` over scored gold spans (mismatch fraction for
start/stop; mean IoU shortfall for `iou`). `unique_rate` / `shared_rate` split
mismatches into tool-only vs same-wrong-coordinate as the other family
(IgBLAST vs SWIG). Per-tool cells still include IoU; joint gallery cells are
start/stop only and carry a pairwise closer matrix.
"""
function span_gallery_cell(pred::AbstractVector{CallRecord},
                           gold::AbstractVector{CallRecord},
                           idx::AbstractDict{String,String},
                           panel::AbstractString,
                           tool::AbstractString,
                           locus,
                           kind;
                           n_sample::Integer = SPAN_SAMPLE_N,
                           rng::AbstractRNG = Random.default_rng(),
                           tool_order = String[String(tool)])
    span_gallery_cell(Dict{String,Vector{CallRecord}}(String(tool) => pred),
                      gold, idx, panel, tool, locus, kind; n_sample, rng, tool_order)
end

function span_gallery_cell(preds::AbstractDict{<:AbstractString,<:AbstractVector{CallRecord}},
                           gold::AbstractVector{CallRecord},
                           idx::AbstractDict{String,String},
                           panel::AbstractString,
                           tool::AbstractString,
                           locus,
                           kind;
                           n_sample::Integer = SPAN_SAMPLE_N,
                           rng::AbstractRNG = Random.default_rng(),
                           tool_order = sort!(collect(String.(keys(preds)))))
    focus = String(tool)
    pred = preds[focus]
    names = String[String(t) for t in tool_order if haskey(preds, t)]
    dstart = Dict{String,Vector{Int}}(t => Int[] for t in names)
    dstop = Dict{String,Vector{Int}}(t => Int[] for t in names)
    miss_start = Dict{String,Int}(t => 0 for t in names)
    miss_stop = Dict{String,Int}(t => 0 for t in names)
    n = 0
    n_bad = 0
    n_shared = 0
    mass = 0.0
    unique_ix = Int[]
    shared_ix = Int[]
    share_of = Dict{Int,String}()
    for i in eachindex(gold)
        gspan = record_span(gold[i], locus)
        isempty(gspan) && continue
        n += 1
        for t in names
            pspan = record_span(preds[t][i], locus)
            ds = span_delta(Val(:start), pspan, gspan)
            dp = span_delta(Val(:stop), pspan, gspan)
            if isnothing(ds)
                miss_start[t] += 1
            else
                push!(dstart[t], ds)
            end
            if isnothing(dp)
                miss_stop[t] += 1
            else
                push!(dstop[t], dp)
            end
        end
        pspan = record_span(pred[i], locus)
        sc = Float64(kind_score(kind, pspan, gspan))
        mass += 1 - sc
        kind_mismatch(kind, pspan, gspan) || continue
        n_bad += 1
        coord = span_coord(kind, pspan)
        shared = other_family_shares(preds, names, focus, i, locus, kind, coord)
        if shared
            n_shared += 1
            push!(shared_ix, i)
            share_of[i] = "shared"
        else
            push!(unique_ix, i)
            share_of[i] = "unique"
        end
    end
    chosen = sample_unique_shared(unique_ix, shared_ix, n_sample, rng)
    samples = Dict{String,Any}[]
    for i in chosen
        peers = Tuple{String,CallRecord}[]
        for t in names
            push!(peers, (t, preds[t][i]))
        end
        push!(samples, sample_record(gold[i], focus, pred[i], peers, locus, idx;
                                     share = get(share_of, i, "")))
    end
    n_unique = n_bad - n_shared
    offsets = Dict{String,Any}()
    for t in names
        offsets[t] = Dict{String,Any}(
            "start" => compact_offset_hist(dstart[t], miss_start[t]),
            "stop" => compact_offset_hist(dstop[t], miss_stop[t]),
        )
    end
    Dict{String,Any}(
        "panel" => String(panel),
        "pred" => focus,
        "ref" => ":gold",
        "locus" => locus_key(locus),
        "kind" => kind_key(kind),
        "n_scored" => n,
        "n_disagree" => n_bad,
        "n_unique" => n_unique,
        "n_shared" => n_shared,
        "n_sample" => length(samples),
        "rate" => n == 0 ? NaN : mass / n,
        "unique_rate" => n == 0 ? NaN : n_unique / n,
        "shared_rate" => n == 0 ? NaN : n_shared / n,
        "offsets" => offsets,
        "samples" => samples,
    )
end

function span_gallery(pred::AbstractVector{CallRecord},
                      gold::AbstractVector{CallRecord},
                      gp::GermlinePaths,
                      panel::AbstractString,
                      tool::AbstractString;
                      n_sample::Integer = SPAN_SAMPLE_N,
                      rng::AbstractRNG = Random.default_rng(),
                      tool_order = String[String(tool)])
    span_gallery(Dict{String,Vector{CallRecord}}(String(tool) => pred),
                 gold, gp, panel, tool; n_sample, rng, tool_order)
end

function span_gallery(preds::AbstractDict{<:AbstractString,<:AbstractVector{CallRecord}},
                      gold::AbstractVector{CallRecord},
                      gp::GermlinePaths,
                      panel::AbstractString,
                      tool::AbstractString;
                      n_sample::Integer = SPAN_SAMPLE_N,
                      rng::AbstractRNG = Random.default_rng(),
                      tool_order = sort!(collect(String.(keys(preds)))))
    idx = germline_sequence_index(gp)
    cells = Dict{String,Any}[]
    for locus in (Val(:v), Val(:d), Val(:j))
        for kind in (Val(:start), Val(:stop), Val(:iou))
            push!(cells, span_gallery_cell(preds, gold, idx, panel, tool, locus, kind;
                                           n_sample, rng, tool_order))
        end
    end
    cells
end

"""One cell per locus × start/stop with every tool's rates, closer matrix, and MSA samples."""
function span_gallery(preds::AbstractDict{<:AbstractString,<:AbstractVector{CallRecord}},
                      gold::AbstractVector{CallRecord},
                      gp::GermlinePaths,
                      panel::AbstractString;
                      n_sample::Integer = SPAN_SAMPLE_N,
                      rng::AbstractRNG = Random.default_rng(),
                      tool_order = sort!(collect(String.(keys(preds)))))
    idx = germline_sequence_index(gp)
    names = String[String(t) for t in tool_order if haskey(preds, t)]
    cells = Dict{String,Any}[]
    for locus in (Val(:v), Val(:d), Val(:j))
        for kind in (Val(:start), Val(:stop))
            push!(cells, span_gallery_joint_cell(preds, gold, idx, panel, names, locus, kind;
                                                 n_sample, rng))
        end
    end
    cells
end

function span_gallery_joint_cell(preds, gold, idx, panel, names, locus, kind;
                                 n_sample::Integer, rng::AbstractRNG)
    nt = length(names)
    n = 0
    dstart = Dict{String,Vector{Int}}(t => Int[] for t in names)
    dstop = Dict{String,Vector{Int}}(t => Int[] for t in names)
    miss_start = Dict{String,Int}(t => 0 for t in names)
    miss_stop = Dict{String,Int}(t => 0 for t in names)
    n_bad = Dict{String,Int}(t => 0 for t in names)
    n_exact_tool = Dict{String,Int}(t => 0 for t in names)
    mass = Dict{String,Float64}(t => 0.0 for t in names)
    win_ix = [Int[] for _ in 1:nt, _ in 1:nt]
    n_win_exact = zeros(Int, nt, nt)
    n_win_closer = zeros(Int, nt, nt)
    n_tie = zeros(Int, nt, nt)
    for i in eachindex(gold)
        gspan = record_span(gold[i], locus)
        isempty(gspan) && continue
        n += 1
        errs = Vector{Int}(undef, nt)
        for (ti, t) in enumerate(names)
            pspan = record_span(preds[t][i], locus)
            ds = span_delta(Val(:start), pspan, gspan)
            dp = span_delta(Val(:stop), pspan, gspan)
            if isnothing(ds)
                miss_start[t] += 1
            else
                push!(dstart[t], ds)
            end
            if isnothing(dp)
                miss_stop[t] += 1
            else
                push!(dstop[t], dp)
            end
            sc = Float64(kind_score(kind, pspan, gspan))
            mass[t] += 1 - sc
            kind_mismatch(kind, pspan, gspan) && (n_bad[t] += 1)
            errs[ti] = edge_abs_error(kind, pspan, gspan)
            errs[ti] == 0 && (n_exact_tool[t] += 1)
        end
        for a in 1:nt
            for b in 1:nt
                a == b && continue
                if errs[a] < errs[b]
                    push!(win_ix[a, b], i)
                    if errs[a] == 0
                        n_win_exact[a, b] += 1
                    else
                        n_win_closer[a, b] += 1
                    end
                elseif errs[a] == errs[b]
                    n_tie[a, b] += 1
                end
            end
        end
    end
    pair_n = Tuple{Int,Int,Int}[(a, b, length(win_ix[a, b])) for a in 1:nt for b in 1:nt if a != b]
    sort!(pair_n; by = p -> (-p[3], names[p[1]], names[p[2]]))
    contrasts = Dict{String,Any}[]
    for (rank, (a, b, nwin)) in enumerate(pair_n)
        push!(contrasts, build_contrast_pack(preds, gold, idx, names, win_ix[a, b],
                                             names[a], names[b],
                                             n_win_exact[a, b], n_win_closer[a, b],
                                             locus, kind, n_sample, rng;
                                             store_samples = contrast_sample_ok(nwin, n, rank)))
    end
    win_n = [a == b ? 0 : length(win_ix[a, b]) for a in 1:nt, b in 1:nt]
    closer = Dict{String,Any}(
        "tools" => names,
        "exact" => Any[n == 0 ? nothing : n_exact_tool[t] / n for t in names],
        "win" => pair_frac(win_n, n),
        "tie" => pair_frac(n_tie, n),
    )
    offsets = Dict{String,Any}()
    tool_stats = Dict{String,Any}()
    for t in names
        offsets[t] = Dict{String,Any}(
            "start" => compact_offset_hist(dstart[t], miss_start[t]),
            "stop" => compact_offset_hist(dstop[t], miss_stop[t]),
        )
        tool_stats[t] = Dict{String,Any}(
            "n_disagree" => n_bad[t],
            "n_exact" => n_exact_tool[t],
            "rate" => n == 0 ? NaN : mass[t] / n,
        )
    end
    Dict{String,Any}(
        "panel" => String(panel),
        "locus" => locus_key(locus),
        "kind" => kind_key(kind),
        "n_scored" => n,
        "n_sample" => Int(n_sample),
        "tools" => tool_stats,
        "offsets" => offsets,
        "closer" => closer,
        "contrasts" => contrasts,
    )
end
