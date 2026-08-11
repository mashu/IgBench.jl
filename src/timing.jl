# timing.jl — Wall-clock throughput around annotate.

"""How to time annotators on a fixed panel."""
struct TimingSpec
    warmup::Int
    repeats::Int
    threads::Int
end

TimingSpec(; warmup::Integer = 1, repeats::Integer = 3, threads::Integer = 1) =
    TimingSpec(Int(warmup), Int(repeats), Int(threads))

"""Per-tool timing summary."""
struct TimingResult
    tool::String
    wall_s::Vector{Float64}
    n_sequences::Int
    seq_per_s::Float64
end

function timing_to_dict(t::TimingResult)
    Dict{String,Any}(
        "tool" => t.tool,
        "wall_s" => t.wall_s,
        "n_sequences" => t.n_sequences,
        "seq_per_s" => t.seq_per_s,
        "median_wall_s" => isempty(t.wall_s) ? NaN : median(t.wall_s),
    )
end

"""
Run warmup + timed repeats of `annotate` on the same inputs.
Returns [`TimingResult`](@ref); `repeats == 0` skips timing.
"""
function time_annotate(tool::AbstractAnnotator, sequences, ids,
                       germline::GermlinePaths, spec::TimingSpec; kwargs...)
    n = length(ids)
    for _ in 1:spec.warmup
        annotate(tool, sequences, ids, germline; kwargs...)
    end
    spec.repeats <= 0 && return TimingResult(tool_name(tool), Float64[], n, NaN)
    walls = Float64[]
    for _ in 1:spec.repeats
        t0 = time_ns()
        annotate(tool, sequences, ids, germline; kwargs...)
        push!(walls, (time_ns() - t0) / 1e9)
    end
    med = median(walls)
    sps = med > 0 ? n / med : Inf
    TimingResult(tool_name(tool), walls, n, sps)
end
