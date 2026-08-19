# timing.jl — Wall-clock throughput of the scored annotate.

"""How to time annotators on a fixed panel."""
struct TimingSpec
    warmup::Int
    repeats::Int
    threads::Int
end

TimingSpec(; warmup::Integer = 0, repeats::Integer = 1, threads::Integer = 1) =
    TimingSpec(Int(warmup), Int(repeats), Int(threads))

"""Per-tool timing of the scored annotate pass."""
struct TimingResult
    tool::String
    wall_s::Float64
    n_sequences::Int
    seq_per_s::Float64
end

function timing_to_dict(t::TimingResult)
    Dict{String,Any}(
        "tool" => t.tool,
        "wall_s" => t.wall_s,
        "n_sequences" => t.n_sequences,
        "seq_per_s" => t.seq_per_s,
    )
end

"""
Run `annotate` once and record wall-clock seconds.

This is the scored pass: callers must use the returned records for metrics.
`spec.repeats <= 0` still annotates but returns `seq_per_s = NaN` (untimed).
`spec.warmup` is ignored (kept on the struct for API stability).
"""
function annotate_timed(tool::AbstractAnnotator, sequences, ids,
                        germline::GermlinePaths, spec::TimingSpec; kwargs...)
    n = length(ids)
    t0 = time_ns()
    rows = annotate(tool, sequences, ids, germline; kwargs...)
    wall = (time_ns() - t0) / 1e9
    if spec.repeats <= 0
        return rows, TimingResult(tool_name(tool), wall, n, NaN)
    end
    sps = wall > 0 ? n / wall : Inf
    rows, TimingResult(tool_name(tool), wall, n, sps)
end

"""Standalone timed annotate (discards records). Prefer [`annotate_timed`](@ref)."""
function time_annotate(tool::AbstractAnnotator, sequences, ids,
                       germline::GermlinePaths, spec::TimingSpec; kwargs...)
    _rows, tr = annotate_timed(tool, sequences, ids, germline, spec; kwargs...)
    tr
end
