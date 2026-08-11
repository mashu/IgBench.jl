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
    println(io, "| panel | tool | seq/s | n |")
    println(io, "|---|---|---:|---:|")
    for t in timing
        sps = t["seq_per_s"]
        sps_s = isnan(Float64(sps)) ? "—" : string(round(Float64(sps); digits = 1))
        println(io, "| ", t["panel"], " | ", t["tool"], " | ", sps_s, " | ",
                t["n_sequences"], " |")
    end
    String(take!(io))
end

write_summary_if_full(::DiagnosticMode, store, name, metrics, timing) = nothing
write_summary_if_full(::FullReportMode, store, name, metrics, timing) =
    write_summary!(store, format_summary(name, metrics, timing))
