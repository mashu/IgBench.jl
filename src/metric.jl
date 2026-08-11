# metric.jl — Pluggable accuracy / agreement metrics.

"""Metric comparing two aligned [`CallRecord`](@ref) vectors."""
abstract type AbstractMetric end

metric_name(m::AbstractMetric) = string(nameof(typeof(m)))

"""Evaluate `pred` against `ref` (same length / order)."""
function evaluate end

function locus_scores(pred::AbstractVector{CallRecord}, ref::AbstractVector{CallRecord},
                      match_fn)
    length(pred) == length(ref) || error("pred/ref length mismatch")
    n = length(pred)
    nv = nd = nj = 0
    dn = 0
    for i in 1:n
        match_fn(pred[i].v_call, ref[i].v_call) && (nv += 1)
        match_fn(pred[i].j_call, ref[i].j_call) && (nj += 1)
        gold_d = !isempty(strip(ref[i].d_call))
        if gold_d
            dn += 1
            match_fn(pred[i].d_call, ref[i].d_call) && (nd += 1)
        end
    end
    MetricValue(v = n == 0 ? NaN : nv / n,
                d = dn == 0 ? NaN : nd / dn,
                j = n == 0 ? NaN : nj / n,
                n = n, d_n = dn)
end

struct ExactCallAccuracy <: AbstractMetric end
metric_name(::ExactCallAccuracy) = "exact"
evaluate(::ExactCallAccuracy, pred, ref) = locus_scores(pred, ref, exact_call_match)

struct AlleleCallAccuracy <: AbstractMetric end
metric_name(::AlleleCallAccuracy) = "allele"
evaluate(::AlleleCallAccuracy, pred, ref) = locus_scores(pred, ref, allele_call_match)

struct GeneCallAccuracy <: AbstractMetric end
metric_name(::GeneCallAccuracy) = "gene"
evaluate(::GeneCallAccuracy, pred, ref) = locus_scores(pred, ref, gene_call_match)

"""Mean span IoU per locus; rows with empty either-side span skipped."""
struct SpanIoU <: AbstractMetric end
metric_name(::SpanIoU) = "span_iou"

function evaluate(::SpanIoU, pred::AbstractVector{CallRecord},
                  ref::AbstractVector{CallRecord})
    length(pred) == length(ref) || error("pred/ref length mismatch")
    function mean_iou(getp, getr)
        vals = Float64[]
        for i in eachindex(pred)
            u = span_iou(getp(pred[i]), getr(ref[i]))
            isnan(u) || push!(vals, u)
        end
        isempty(vals) ? NaN : sum(vals) / length(vals)
    end
    MetricValue(
        v = mean_iou(r -> r.v_span, r -> r.v_span),
        d = mean_iou(r -> r.d_span, r -> r.d_span),
        j = mean_iou(r -> r.j_span, r -> r.j_span),
        n = length(pred),
        d_n = count(r -> !isempty(r.d_span), ref),
    )
end

"""Default metric set for tool vs gold / tool vs tool."""
default_metrics() = AbstractMetric[ExactCallAccuracy(), AlleleCallAccuracy(),
                                   GeneCallAccuracy(), SpanIoU()]

agreement_metrics() = AbstractMetric[AlleleCallAccuracy(), GeneCallAccuracy()]
