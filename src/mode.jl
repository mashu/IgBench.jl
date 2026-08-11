# mode.jl — Diagnostic vs full-report execution policy.

abstract type RunMode end

"""Cheap mid-training eval: small N, light/no timing, predictions off by default."""
struct DiagnosticMode <: RunMode
    max_sequences::Int
    store_predictions::Bool
    timing_repeats::Int
end

DiagnosticMode(; max_sequences::Integer = 512,
               store_predictions::Bool = false,
               timing_repeats::Integer = 0) =
    DiagnosticMode(Int(max_sequences), store_predictions, Int(timing_repeats))

"""End-of-train / standalone report: full N, timings, predictions on."""
struct FullReportMode <: RunMode
    store_predictions::Bool
    timing_warmup::Int
    timing_repeats::Int
end

FullReportMode(; store_predictions::Bool = true,
               timing_warmup::Integer = 1,
               timing_repeats::Integer = 3) =
    FullReportMode(store_predictions, Int(timing_warmup), Int(timing_repeats))

mode_name(::DiagnosticMode) = "diagnostic"
mode_name(::FullReportMode) = "full"

store_predictions(m::DiagnosticMode) = m.store_predictions
store_predictions(m::FullReportMode) = m.store_predictions

function timing_spec(m::DiagnosticMode)
    TimingSpec(warmup = 0, repeats = m.timing_repeats, threads = 1)
end

function timing_spec(m::FullReportMode)
    TimingSpec(warmup = m.timing_warmup, repeats = m.timing_repeats, threads = 1)
end

"""Cap panel size under diagnostic mode; full mode keeps requested n."""
effective_n(m::DiagnosticMode, requested::Integer) = min(Int(requested), m.max_sequences)
effective_n(::FullReportMode, requested::Integer) = Int(requested)

merge_timing(m::DiagnosticMode, base::TimingSpec) =
    TimingSpec(warmup = 0, repeats = m.timing_repeats, threads = base.threads)

merge_timing(m::FullReportMode, base::TimingSpec) =
    TimingSpec(warmup = m.timing_warmup, repeats = m.timing_repeats, threads = base.threads)
