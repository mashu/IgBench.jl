# call_events.jl — Gold span lengths for missed (gold call, empty pred) reads.

function compact_count_hist(values::AbstractVector{Int})
    acc = Dict{Int,Int}()
    for x in values
        acc[x] = get(acc, x, 0) + 1
    end
    ks = sort!(collect(keys(acc)))
    bins = Vector{Vector{Int}}(undef, length(ks))
    for i in eachindex(ks)
        k = ks[i]
        bins[i] = Int[k, acc[k]]
    end
    Dict{String,Any}("n" => length(values), "bins" => bins)
end

function miss_read_record(gold::CallRecord, sp::Span, len::Int)
    seq = gold.sequence
    Dict{String,Any}(
        "sequence_id" => gold.sequence_id,
        "sequence" => seq,
        "length" => len,
        "v_call" => gold.v_call,
        "d_call" => gold.d_call,
        "j_call" => gold.j_call,
        "v_nt" => span_nt(seq, gold.v_span),
        "d_nt" => span_nt(seq, gold.d_span),
        "j_nt" => span_nt(seq, gold.j_span),
        "span" => isempty(sp) ? Int[] : Int[sp.start, sp.stop],
    )
end

function collect_miss_reads(pred, gold, pred_call, gold_call, gold_span)
    n = 0
    lens = Int[]
    reads = Dict{String,Any}[]
    for i in eachindex(gold)
        call_field_empty(gold_call(gold[i])) && continue
        call_field_empty(pred_call(pred[i])) || continue
        n += 1
        g = gold[i]
        sp = gold_span(g)
        len = isempty(sp) ? 0 : length(sp)
        isempty(sp) || push!(lens, len)
        push!(reads, miss_read_record(g, sp, len))
    end
    n, lens, reads
end

"""Gold span-length histograms plus FASTA-ready reads for missed calls."""
function call_miss_hists(pred::AbstractVector{CallRecord}, gold::AbstractVector{CallRecord},
                         panel, tool)
    length(pred) == length(gold) || error("pred/ref length mismatch")
    nv, vlen, vreads = collect_miss_reads(pred, gold, r -> r.v_call, r -> r.v_call, r -> r.v_span)
    nd, dlen, dreads = collect_miss_reads(pred, gold, r -> r.d_call, r -> r.d_call, r -> r.d_span)
    nj, jlen, jreads = collect_miss_reads(pred, gold, r -> r.j_call, r -> r.j_call, r -> r.j_span)
    Dict{String,Any}(
        "panel" => String(panel),
        "pred" => String(tool),
        "ref" => ":gold",
        "v_miss" => nv,
        "d_miss" => nd,
        "j_miss" => nj,
        "v" => compact_count_hist(vlen),
        "d" => compact_count_hist(dlen),
        "j" => compact_count_hist(jlen),
        "v_reads" => vreads,
        "d_reads" => dreads,
        "j_reads" => jreads,
    )
end
