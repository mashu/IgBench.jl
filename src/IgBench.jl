"""
    IgBench

Extensible VDJ benchmark harness: pluggable annotators, metrics, and panels.
Library-first API ([`run_suite`](@ref)); optional standalone `scripts/run_bench.jl`.
IgBLAST.jl is the only built-in tool package dependency — SwiftIG / IgFormer plug in
via [`SwiftIGAnnotator`](@ref) / [`CallableAnnotator`](@ref).
"""
module IgBench

using Random
using Statistics
using Dates
using JSON
using CodecZlib
using IgSim

include("types.jl")
include("match.jl")
include("airr_io.jl")
include("manifest.jl")
include("annotator.jl")
include("igblast.jl")
include("swiftig.jl")
include("metric.jl")
include("timing.jl")
include("mode.jl")
include("panel.jl")
include("store.jl")
include("report.jl")
include("report_html.jl")
include("suite.jl")
include("defaults.jl")

export
    Span, EMPTY_SPAN, GermlinePaths, CallRecord, MetricValue,
    AirrSource, SimSource, DatasetManifest, manifest_to_dict,
    AbstractAnnotator, tool_name, annotate,
    CallableAnnotator, FakeAnnotator,
    IgBLASTAnnotator, SwiftIGAnnotator,
    AbstractMetric, metric_name, evaluate,
    AlleleAccuracy, SpanIoU, SpanExact, SpanStart, SpanStop,
    default_metrics, agreement_metrics, call_metrics, call_records,
    TimingSpec, TimingResult, time_annotate, annotate_timed,
    RunMode, DiagnosticMode, FullReportMode, mode_name,
    AbstractPanel, PanelData, SimGoldPanel, AirrPanel, load_panel,
    PanelCache, load_panel_cached, SIM_SETS,
    AbstractRunStore, DirectoryRunStore, NullRunStore, SCHEMA_VERSION,
    CompareSpec, BenchSuite, BenchResult, run_suite,
    default_compares, suite_from_manifest,
    read_airr_calls, write_airr_calls,
    parse_allele_calls, normalize_allele, allele_gene, call_field_empty,
    primary_allele_call, allele_score,
    span_iou, span_exact, span_start, span_stop,
    germline_allele_names, validate_gold_calls!

end # module
