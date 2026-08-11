# match.jl — Allele / gene call matching (AIRR multi-call aware).

"""Split IgBLAST-style multi-calls on `/` and `,`."""
function parse_allele_calls(field::AbstractString)
    s = strip(String(field))
    (isempty(s) || s == "NA") && return String[]
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

"""Exact string equality (no splitting / normalization)."""
exact_call_match(pred::AbstractString, gold::AbstractString) =
    String(pred) == String(gold)

"""True if any pred token matches any gold token at allele level."""
function allele_call_match(pred::AbstractString, gold::AbstractString)
    golds = parse_allele_calls(gold)
    preds = parse_allele_calls(pred)
    isempty(golds) && return isempty(preds)
    isempty(preds) && return false
    gset = Set(normalize_allele(g) for g in golds)
    any(normalize_allele(p) in gset for p in preds)
end

"""True if any pred token matches any gold token at gene level."""
function gene_call_match(pred::AbstractString, gold::AbstractString)
    golds = parse_allele_calls(gold)
    preds = parse_allele_calls(pred)
    isempty(golds) && return isempty(preds)
    isempty(preds) && return false
    gset = Set(allele_gene(g) for g in golds)
    any(allele_gene(p) in gset for p in preds)
end

"""Intersection-over-union of two spans; `NaN` if either empty."""
function span_iou(a::Span, b::Span)
    (isempty(a) || isempty(b)) && return NaN
    inter = max(0, min(a.stop, b.stop) - max(a.start, b.start) + 1)
    union = length(a) + length(b) - inter
    union == 0 ? NaN : inter / union
end
