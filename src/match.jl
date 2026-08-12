# match.jl — Allele / gene call matching (AIRR multi-call aware).
#
# Canonical tool-vs-gold scoring rules (see also [`FractionalCallAccuracy`](@ref)):
# 1. Empty / missing / `NA` / `.` **gold** → skip (not in the denominator).
# 2. Non-empty gold + empty **pred** → score `0` (false negative).
# 3. Multi-allelic pred: any token matching gold → `1/n` (`FractionalCallAccuracy`);
#    full-string match → `1` (preserves IMGT dual names containing `/`).

"""True for blank, `NA`, or `.` call fields."""
function call_field_empty(field::AbstractString)
    s = strip(String(field))
    isempty(s) || s == "NA" || s == "."
end

call_field_empty(::Missing) = true

"""Split IgBLAST-style multi-calls on `/` and `,`."""
function parse_allele_calls(field::AbstractString)
    s = strip(String(field))
    (isempty(s) || s == "NA" || s == ".") && return String[]
    out = String[]
    for p in split(s, r"[/,]")
        t = String(strip(p))
        isempty(t) || push!(out, t)
    end
    unique!(out)
    out
end

"""Strip personalized `_S####` suffixes."""
normalize_allele(call::AbstractString) = replace(String(strip(call)), r"_S\d+$" => "")

"""Gene locus before `*` (after `_S*` strip)."""
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

"""Exact string equality (no splitting / normalization)."""
exact_call_match(pred::AbstractString, gold::AbstractString) =
    String(pred) == String(gold)

exact_call_match(::Missing, gold::AbstractString) = exact_call_match("", gold)
exact_call_match(pred::AbstractString, ::Missing) = exact_call_match(pred, "")
exact_call_match(::Missing, ::Missing) = true

"""True if any pred token matches any gold token at allele level (full point)."""
function allele_call_match(pred::AbstractString, gold::AbstractString)
    golds = parse_allele_calls(gold)
    preds = parse_allele_calls(pred)
    isempty(golds) && return isempty(preds)
    isempty(preds) && return false
    gset = Set(normalize_allele(g) for g in golds)
    any(normalize_allele(p) in gset for p in preds)
end

allele_call_match(::Missing, gold::AbstractString) = allele_call_match("", gold)
allele_call_match(pred::AbstractString, ::Missing) = allele_call_match(pred, "")
allele_call_match(::Missing, ::Missing) = true

"""True if any pred token matches any gold token at gene level."""
function gene_call_match(pred::AbstractString, gold::AbstractString)
    golds = parse_allele_calls(gold)
    preds = parse_allele_calls(pred)
    isempty(golds) && return isempty(preds)
    isempty(preds) && return false
    gset = Set(allele_gene(g) for g in golds)
    any(allele_gene(p) in gset for p in preds)
end

gene_call_match(::Missing, gold::AbstractString) = gene_call_match("", gold)
gene_call_match(pred::AbstractString, ::Missing) = gene_call_match(pred, "")
gene_call_match(::Missing, ::Missing) = true

"""
Fractional multi-call score in `[0, 1]`.

Empty gold → `0` (caller should skip). Empty pred → `0`. Full-string match → `1`.
Else if any `,`/`/`-split pred token matches any gold token → `1/n`.
"""
function fractional_call_score(pred::AbstractString, gold::AbstractString)
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

fractional_call_score(pred::Missing, gold::AbstractString) =
    fractional_call_score("", gold)
fractional_call_score(pred::AbstractString, gold::Missing) =
    fractional_call_score(pred, "")
fractional_call_score(::Missing, ::Missing) = 0.0

"""First-comma pred equals full gold string (ablation / legacy)."""
function primary_call_match(pred::AbstractString, gold::AbstractString)
    call_field_empty(gold) && return false
    primary_allele_call(pred) == strip(String(gold))
end

"""Intersection-over-union of two spans; `NaN` if either empty."""
function span_iou(a::Span, b::Span)
    (isempty(a) || isempty(b)) && return NaN
    inter = max(0, min(a.stop, b.stop) - max(a.start, b.start) + 1)
    union = length(a) + length(b) - inter
    union == 0 ? NaN : inter / union
end
