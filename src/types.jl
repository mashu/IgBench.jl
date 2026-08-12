# types.jl — Core records shared by panels, tools, and metrics.

"""1-based inclusive sequence interval; empty when `stop < start`."""
struct Span
    start::Int
    stop::Int
end

const EMPTY_SPAN = Span(0, -1)
Base.isempty(s::Span) = s.stop < s.start
Base.length(s::Span) = isempty(s) ? 0 : s.stop - s.start + 1

"""Paths to V/D/J germline FASTA files (`d` may be `nothing`)."""
struct GermlinePaths
    v::String
    d::Union{String,Nothing}
    j::String
end

GermlinePaths(; v::AbstractString, d = nothing, j::AbstractString) =
    GermlinePaths(String(v), isnothing(d) ? nothing : String(d), String(j))

"""
One annotated read: sequence plus V/D/J calls and optional coordinate spans.
"""
struct CallRecord
    sequence_id::String
    sequence::String
    v_call::String
    d_call::String
    j_call::String
    v_span::Span
    d_span::Span
    j_span::Span
end

CallRecord(id, seq, v, d, j) =
    CallRecord(String(id), String(seq), String(v), String(d), String(j),
               EMPTY_SPAN, EMPTY_SPAN, EMPTY_SPAN)

"""Scalar bundle returned by metrics (JSON-serializable)."""
struct MetricValue
    v::Float64
    d::Float64
    j::Float64
    n::Int
    d_n::Int
    v_n::Int
    j_n::Int
end

MetricValue(; v, d, j, n, d_n = n, v_n = n, j_n = n) =
    MetricValue(Float64(v), Float64(d), Float64(j), Int(n), Int(d_n), Int(v_n), Int(j_n))

function metric_dict(m::MetricValue)
    Dict{String,Any}("v" => m.v, "d" => m.d, "j" => m.j, "n" => m.n,
                     "d_n" => m.d_n, "v_n" => m.v_n, "j_n" => m.j_n)
end
