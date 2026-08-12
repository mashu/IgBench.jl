# metric.jl — Pluggable accuracy / agreement metrics.
#
# Prefer [`FractionalCallAccuracy`](@ref) for tool-vs-gold benches: fair `1/n`
# credit on multi-calls, empty-gold rows skipped, empty pred = false negative.

"""Metric comparing two aligned [`CallRecord`](@ref) vectors."""
abstract type AbstractMetric end

metric_name(m::AbstractMetric) = string(nameof(typeof(m)))

"""Evaluate `pred` against `ref` (same length / order)."""
function evaluate end

"""
Per-locus mean of `score_fn(pred_call, gold_call)` over non-empty gold.

`score_fn` may return `Bool` or `Real` in `[0, 1]`. Empty gold is skipped on
V, D, and J. Empty pred vs present gold should score `0` / `false`.
"""
function locus_scores(pred::AbstractVector{CallRecord}, ref::AbstractVector{CallRecord},
                      score_fn)
    length(pred) == length(ref) || error("pred/ref length mismatch")
    n = length(pred)
    sv = sd = sj = 0.0
    vn = dn = jn = 0
    for i in 1:n
        if !call_field_empty(ref[i].v_call)
            vn += 1
            sv += Float64(score_fn(pred[i].v_call, ref[i].v_call))
        end
        if !call_field_empty(ref[i].j_call)
            jn += 1
            sj += Float64(score_fn(pred[i].j_call, ref[i].j_call))
        end
        if !call_field_empty(ref[i].d_call)
            dn += 1
            sd += Float64(score_fn(pred[i].d_call, ref[i].d_call))
        end
    end
    MetricValue(v = vn == 0 ? NaN : sv / vn,
                d = dn == 0 ? NaN : sd / dn,
                j = jn == 0 ? NaN : sj / jn,
                n = n, d_n = dn, v_n = vn, j_n = jn)
end

"""Build aligned [`CallRecord`](@ref)s from parallel V/D/J call columns."""
function call_records(v_calls::AbstractVector, d_calls::AbstractVector,
                      j_calls::AbstractVector; prefix::AbstractString = "r")
    n = length(v_calls)
    length(d_calls) == n && length(j_calls) == n ||
        error("V/D/J call column length mismatch")
    CallRecord[CallRecord("$prefix$i", "", string(v_calls[i]), string(d_calls[i]),
                          string(j_calls[i])) for i in 1:n]
end

"""
**Standard** tool-vs-gold metric: fractional `1/n` multi-call credit.

See [`fractional_call_score`](@ref). Empty gold skipped; empty pred = 0.
"""
struct FractionalCallAccuracy <: AbstractMetric end
metric_name(::FractionalCallAccuracy) = "fractional"
evaluate(::FractionalCallAccuracy, pred, ref) =
    locus_scores(pred, ref, fractional_call_score)

"""Full-string exact match (after empty-gold skip)."""
struct ExactCallAccuracy <: AbstractMetric end
metric_name(::ExactCallAccuracy) = "exact"
evaluate(::ExactCallAccuracy, pred, ref) = locus_scores(pred, ref, exact_call_match)

"""Any-token allele match → full point (lenient multi-call; not `1/n`)."""
struct AlleleCallAccuracy <: AbstractMetric end
metric_name(::AlleleCallAccuracy) = "allele"
evaluate(::AlleleCallAccuracy, pred, ref) = locus_scores(pred, ref, allele_call_match)

"""Gene-level match."""
struct GeneCallAccuracy <: AbstractMetric end
metric_name(::GeneCallAccuracy) = "gene"
evaluate(::GeneCallAccuracy, pred, ref) = locus_scores(pred, ref, gene_call_match)

"""First comma-separated call vs full gold (ablation)."""
struct PrimaryCallAccuracy <: AbstractMetric end
metric_name(::PrimaryCallAccuracy) = "primary"
evaluate(::PrimaryCallAccuracy, pred, ref) = locus_scores(pred, ref, primary_call_match)

"""Mean span IoU per locus; rows with empty either-side span skipped."""
struct SpanIoU <: AbstractMetric end
metric_name(::SpanIoU) = "span_iou"

function evaluate(::SpanIoU, pred::AbstractVector{CallRecord},
                  ref::AbstractVector{CallRecord})
    length(pred) == length(ref) || error("pred/ref length mismatch")
    function mean_iou(getp, getr)
        s = 0.0
        n = 0
        for i in eachindex(pred)
            u = span_iou(getp(pred[i]), getr(ref[i]))
            isnan(u) && continue
            s += u
            n += 1
        end
        n == 0 ? NaN : s / n
    end
    MetricValue(
        v = mean_iou(r -> r.v_span, r -> r.v_span),
        d = mean_iou(r -> r.d_span, r -> r.d_span),
        j = mean_iou(r -> r.j_span, r -> r.j_span),
        n = length(pred),
        d_n = count(r -> !isempty(r.d_span), ref),
        v_n = count(r -> !isempty(r.v_span), ref),
        j_n = count(r -> !isempty(r.j_span), ref),
    )
end

"""
Default metric set for tool vs gold / tool vs tool.

Leads with [`FractionalCallAccuracy`](@ref) (canonical bench score).
"""
default_metrics() = AbstractMetric[FractionalCallAccuracy(), ExactCallAccuracy(),
                                   AlleleCallAccuracy(), GeneCallAccuracy(),
                                   PrimaryCallAccuracy(), SpanIoU()]

agreement_metrics() = AbstractMetric[FractionalCallAccuracy(), AlleleCallAccuracy(),
                                     GeneCallAccuracy()]

"""
    call_metrics(v_pred, d_pred, j_pred, v_gold, d_gold, j_gold) -> MetricValue

Convenience wrapper: string columns → [`FractionalCallAccuracy`](@ref).
"""
function call_metrics(v_pred::AbstractVector, d_pred::AbstractVector,
                      j_pred::AbstractVector,
                      v_gold::AbstractVector, d_gold::AbstractVector,
                      j_gold::AbstractVector)
    length(v_pred) == length(v_gold) || error("pred/gold length mismatch")
    evaluate(FractionalCallAccuracy(),
             call_records(v_pred, d_pred, j_pred; prefix = "p"),
             call_records(v_gold, d_gold, j_gold; prefix = "g"))
end
