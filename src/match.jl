# match.jl — Allele scoring and span geometry (AIRR multi-call aware).
#
# Canonical tool-vs-gold allele score ([`AlleleAccuracy`](@ref) / [`allele_score`](@ref)):
# 1. Empty / missing / `NA` / `.` **gold** → skip (not in the denominator).
# 2. Non-empty gold + empty **pred** → score `0` (false negative).
# 3. Full-string equality → `1` (identical multi-calls recapitulate).
# 4. Else split on commas only (slash is part of an IMGT dual name);
#    any pred token equal to any gold token → `1 / n_pred`.

"""True for blank, `NA`, or `.` call fields."""
function call_field_empty(field::AbstractString)
    s = strip(String(field))
    isempty(s) || s == "NA" || s == "."
end

call_field_empty(::Missing) = true

"""
Split IgBLAST-style multi-calls on commas only.

Slash is part of one allele identifier (e.g. `IGHV3-23*01/IGHV3-23D*01`).
"""
function parse_allele_calls(field::AbstractString)
    s = strip(String(field))
    (isempty(s) || s == "NA" || s == ".") && return String[]
    out = String[]
    for p in split(s, ',')
        t = String(strip(p))
        isempty(t) || push!(out, t)
    end
    unique!(out)
    out
end

parse_allele_calls(::Missing) = String[]

"""Strip personalized `_S####` suffixes (analysis helper; not used in allele score)."""
normalize_allele(call::AbstractString) = replace(String(strip(call)), r"_S\d+$" => "")

"""Gene locus before `*` (after `_S*` strip). Analysis helper; not a bench metric."""
function allele_gene(call::AbstractString)
    s = normalize_allele(call)
    i = findfirst('*', s)
    isnothing(i) ? s : s[1:prevind(s, i)]
end

"""
First AIRR allele call (comma-separated list).

Splits on `,` only so IMGT dual names that contain `/` stay intact.
"""
function primary_allele_call(field::AbstractString)
    call_field_empty(field) && return ""
    s = strip(String(field))
    c = findfirst(',', s)
    c === nothing && return s
    String(strip(SubString(s, 1, prevind(s, c))))
end

primary_allele_call(::Missing) = ""

"""
Allele score in `[0, 1]`.

Empty gold → `0` (caller should skip). Empty pred → `0`. Full-string match → `1`.
Else if any comma-split pred token equals any gold token → `1 / n_pred`.
Does not strip `_S` suffixes and does not split on `/`.
"""
function allele_score(pred::AbstractString, gold::AbstractString)
    call_field_empty(gold) && return 0.0
    call_field_empty(pred) && return 0.0
    p = strip(String(pred))
    g = strip(String(gold))
    p == g && return 1.0
    preds = parse_allele_calls(p)
    isempty(preds) && return 0.0
    golds = parse_allele_calls(g)
    isempty(golds) && return 0.0
    gset = Set(golds)
    any(t -> t in gset, preds) || return 0.0
    1.0 / length(preds)
end

allele_score(pred::Missing, gold::AbstractString) = allele_score("", gold)
allele_score(pred::AbstractString, gold::Missing) = allele_score(pred, "")
allele_score(::Missing, ::Missing) = 0.0

"""
Intersection-over-union of pred vs gold (1-based inclusive).

Empty **gold** → `NaN` (caller skips). Empty **pred** vs present gold → `0`.
"""
function span_iou(pred::Span, gold::Span)
    isempty(gold) && return NaN
    isempty(pred) && return 0.0
    inter = max(0, min(pred.stop, gold.stop) - max(pred.start, gold.start) + 1)
    union = length(pred) + length(gold) - inter
    union == 0 ? 0.0 : inter / union
end

"""1 iff pred start and stop both equal gold. Empty gold → `NaN`; empty pred → `0`."""
function span_exact(pred::Span, gold::Span)
    isempty(gold) && return NaN
    isempty(pred) && return 0.0
    (pred.start == gold.start && pred.stop == gold.stop) ? 1.0 : 0.0
end

"""1 iff pred start equals gold. Empty gold → `NaN`; empty pred → `0`."""
function span_start(pred::Span, gold::Span)
    isempty(gold) && return NaN
    isempty(pred) && return 0.0
    pred.start == gold.start ? 1.0 : 0.0
end

"""1 iff pred stop equals gold. Empty gold → `NaN`; empty pred → `0`."""
function span_stop(pred::Span, gold::Span)
    isempty(gold) && return NaN
    isempty(pred) && return 0.0
    pred.stop == gold.stop ? 1.0 : 0.0
end
