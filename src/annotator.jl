# annotator.jl — Pluggable VDJ annotation tools.

"""Tool that maps sequences → [`CallRecord`](@ref) rows."""
abstract type AbstractAnnotator end

"""Stable short name used in compares and artifact filenames."""
function tool_name end

"""
    annotate(tool, sequences, ids, germline; kwargs...) -> Vector{CallRecord}

Must return one record per input id, same order.
"""
function annotate end

"""Wrap any callable `(seqs, ids, germline) -> Vector{CallRecord}`."""
struct CallableAnnotator{F} <: AbstractAnnotator
    name::String
    f::F
end

CallableAnnotator(name::AbstractString, f) = CallableAnnotator(String(name), f)
# support `CallableAnnotator("name") do ... end`
CallableAnnotator(f::Function, name::AbstractString) = CallableAnnotator(String(name), f)

tool_name(a::CallableAnnotator) = a.name

annotate(a::CallableAnnotator, sequences, ids, germline::GermlinePaths; kwargs...) =
    a.f(sequences, ids, germline)

"""
Deterministic stub for tests: copies gold when provided via kwargs, else empty calls.
"""
struct FakeAnnotator <: AbstractAnnotator
    name::String
    gold::Union{Nothing,Vector{CallRecord}}
end

FakeAnnotator(name::AbstractString = "fake"; gold = nothing) =
    FakeAnnotator(String(name), gold)

tool_name(a::FakeAnnotator) = a.name

function annotate(a::FakeAnnotator, sequences, ids, germline::GermlinePaths; kwargs...)
    n = length(ids)
    length(sequences) == n || error("sequences/ids length mismatch")
    if !isnothing(a.gold)
        length(a.gold) == n || error("FakeAnnotator gold length mismatch")
        return a.gold
    end
    CallRecord[CallRecord(ids[i], sequences[i], "", "", "") for i in 1:n]
end
