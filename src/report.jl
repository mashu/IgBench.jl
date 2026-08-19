# report.jl — Human-readable summary.md for full reports.

function format_summary(suite_name::AbstractString,
                        metrics::AbstractVector{<:AbstractDict},
                        timing::AbstractVector{<:AbstractDict})
    io = IOBuffer()
    println(io, "# IgBench report: ", suite_name)
    println(io)
    println(io, "## Metrics")
    println(io)
    println(io, "| panel | pred | ref | metric | V | D | J | n |")
    println(io, "|---|---|---|---:|---:|---:|---:|---:|")
    for r in metrics
        println(io, "| ", r["panel"], " | ", r["pred"], " | ", r["ref"], " | ",
                r["metric"], " | ", round(Float64(r["v"]); digits = 4), " | ",
                round(Float64(r["d"]); digits = 4), " | ",
                round(Float64(r["j"]); digits = 4), " | ", r["n"], " |")
    end
    println(io)
    println(io, "## Timing")
    println(io)
    println(io, "| panel | tool | wall_s | n | seq/s |")
    println(io, "|---|---|---:|---:|---:|")
    for t in timing
        wall = Float64(t["wall_s"])
        sps = t["seq_per_s"]
        wall_s = isnan(wall) ? "—" : string(round(wall; digits = 3))
        sps_s = isnan(Float64(sps)) ? "—" : string(round(Float64(sps); digits = 1))
        println(io, "| ", t["panel"], " | ", t["tool"], " | ", wall_s, " | ",
                t["n_sequences"], " | ", sps_s, " |")
    end
    String(take!(io))
end

write_summary_if_full(::DiagnosticMode, store, name, metrics, timing) = nothing
write_summary_if_full(::FullReportMode, store, name, metrics, timing) =
    write_summary!(store, format_summary(name, metrics, timing))
