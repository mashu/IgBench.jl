# metric.jl — Pluggable accuracy / agreement metrics.
#
# Published tool-vs-gold scores: [`AlleleAccuracy`](@ref), [`CallPresent`](@ref),
# and the span metrics. Empty gold skipped; empty pred vs present gold = 0.

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

"""Like [`locus_scores`](@ref) but only over **empty** gold (spurious pred calls)."""
function locus_absent_gold_scores(pred::AbstractVector{CallRecord},
                                  ref::AbstractVector{CallRecord}, score_fn)
    length(pred) == length(ref) || error("pred/ref length mismatch")
    n = length(pred)
    sv = sd = sj = 0.0
    vn = dn = jn = 0
    for i in 1:n
        if call_field_empty(ref[i].v_call)
            vn += 1
            sv += Float64(score_fn(pred[i].v_call, ref[i].v_call))
        end
        if call_field_empty(ref[i].j_call)
            jn += 1
            sj += Float64(score_fn(pred[i].j_call, ref[i].j_call))
        end
        if call_field_empty(ref[i].d_call)
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
**Standard** tool-vs-gold allele metric: `1/n` on comma-separated ties.

See [`allele_score`](@ref). Empty vs empty = correct; empty gold is not skipped.
"""
struct AlleleAccuracy <: AbstractMetric end
metric_name(::AlleleAccuracy) = "allele"
evaluate(::AlleleAccuracy, pred, ref) = locus_all_scores(pred, ref, allele_score)

"""Per-locus mean of `score_fn` over **all** reads (empty gold is not skipped)."""
function locus_all_scores(pred::AbstractVector{CallRecord}, ref::AbstractVector{CallRecord},
                         score_fn)
    length(pred) == length(ref) || error("pred/ref length mismatch")
    n = length(pred)
    sv = sd = sj = 0.0
    for i in 1:n
        sv += Float64(score_fn(pred[i].v_call, ref[i].v_call))
        sj += Float64(score_fn(pred[i].j_call, ref[i].j_call))
        sd += Float64(score_fn(pred[i].d_call, ref[i].d_call))
    end
    MetricValue(v = n == 0 ? NaN : sv / n,
                d = n == 0 ? NaN : sd / n,
                j = n == 0 ? NaN : sj / n,
                n = n, d_n = n, v_n = n, j_n = n)
end

"""1 if the pred call is non-empty (gold already present). Empty pred = miss."""
call_is_present(pred::AbstractString, ::AbstractString) = !call_field_empty(pred)

"""Fraction of present-gold reads where the tool emitted a call."""
struct CallPresent <: AbstractMetric end
metric_name(::CallPresent) = "call_present"
evaluate(::CallPresent, pred, ref) = locus_scores(pred, ref, call_is_present)

"""Fraction of empty-gold reads where the tool still emitted a call."""
struct CallExtra <: AbstractMetric end
metric_name(::CallExtra) = "call_extra"
evaluate(::CallExtra, pred, ref) = locus_absent_gold_scores(pred, ref, call_is_present)

"""Mean of a span score_fn per locus; empty gold spans skipped."""
function locus_span_scores(pred::AbstractVector{CallRecord},
                           ref::AbstractVector{CallRecord},
                           score_fn)
    length(pred) == length(ref) || error("pred/ref length mismatch")
    function mean_score(getp, getr)
        s = 0.0
        n = 0
        for i in eachindex(pred)
            gold = getr(ref[i])
            isempty(gold) && continue
            s += Float64(score_fn(getp(pred[i]), gold))
            n += 1
        end
        n == 0 ? NaN : s / n
    end
    MetricValue(
        v = mean_score(r -> r.v_span, r -> r.v_span),
        d = mean_score(r -> r.d_span, r -> r.d_span),
        j = mean_score(r -> r.j_span, r -> r.j_span),
        n = length(pred),
        d_n = count(r -> !isempty(r.d_span), ref),
        v_n = count(r -> !isempty(r.v_span), ref),
        j_n = count(r -> !isempty(r.j_span), ref),
    )
end

"""Mean span IoU per locus. Missing pred vs present gold = 0."""
struct SpanIoU <: AbstractMetric end
metric_name(::SpanIoU) = "span_iou"
evaluate(::SpanIoU, pred, ref) = locus_span_scores(pred, ref, span_iou)

"""Fraction of reads where start and stop both equal gold."""
struct SpanExact <: AbstractMetric end
metric_name(::SpanExact) = "span_exact"
evaluate(::SpanExact, pred, ref) = locus_span_scores(pred, ref, span_exact)

"""Fraction of reads where start equals gold."""
struct SpanStart <: AbstractMetric end
metric_name(::SpanStart) = "span_start"
evaluate(::SpanStart, pred, ref) = locus_span_scores(pred, ref, span_start)

"""Fraction of reads where stop equals gold."""
struct SpanStop <: AbstractMetric end
metric_name(::SpanStop) = "span_stop"
evaluate(::SpanStop, pred, ref) = locus_span_scores(pred, ref, span_stop)

"""1 if the pred span is non-empty (gold already present). Empty pred = miss."""
span_is_present(pred::Span, ::Span) = !isempty(pred)

"""Fraction of present-gold spans where the tool emitted an interval."""
struct SpanPresent <: AbstractMetric end
metric_name(::SpanPresent) = "span_present"
evaluate(::SpanPresent, pred, ref) = locus_span_scores(pred, ref, span_is_present)

"""Default metric set: allele score plus independent span metrics."""
default_metrics() = AbstractMetric[AlleleAccuracy(), CallPresent(), CallExtra(), SpanIoU(),
                                   SpanExact(), SpanStart(), SpanStop(), SpanPresent()]

"""Tool-vs-tool: same metric set (allele + spans)."""
agreement_metrics() = default_metrics()

"""
    call_metrics(v_pred, d_pred, j_pred, v_gold, d_gold, j_gold) -> MetricValue

Convenience wrapper: string columns → [`AlleleAccuracy`](@ref).
"""
function call_metrics(v_pred::AbstractVector, d_pred::AbstractVector,
                      j_pred::AbstractVector,
                      v_gold::AbstractVector, d_gold::AbstractVector,
                      j_gold::AbstractVector)
    length(v_pred) == length(v_gold) || error("pred/gold length mismatch")
    evaluate(AlleleAccuracy(),
             call_records(v_pred, d_pred, j_pred; prefix = "p"),
             call_records(v_gold, d_gold, j_gold; prefix = "g"))
end
