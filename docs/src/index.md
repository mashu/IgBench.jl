# IgBench.jl

Extensible V(D)J benchmark harness. Use as a **library** from IgFormer (or
any client) via [`run_suite`](@ref), or standalone via `scripts/run_bench.jl`.

IgBLAST.jl is the only built-in tool package dependency. SwiftIG and IgFormer
register through [`SwiftIGAnnotator`](@ref) / [`CallableAnnotator`](@ref).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgBench.jl")
```

## Library (embedded)

```julia
using IgBench

suite = suite_from_manifest(manifest, tools)
result = run_suite(suite;
                   mode = DiagnosticMode(max_sequences = 256),
                   store = DirectoryRunStore("runs/diagnostic/1000"),
                   step = 1000, tags = (; epoch = 3))
```

## Published metrics

[`AlleleAccuracy`](@ref) (metric name `"allele"`) is the only call score:

1. Empty / `NA` / `.` gold → skip (V, D, and J).
2. Present gold + empty pred → `0`.
3. Full-string match → `1`. Else comma-split (slash is part of an IMGT dual
   name); any matching token → `1 / n_pred`.

Segmentation (independent of allele): [`SpanIoU`](@ref), [`SpanExact`](@ref),
[`SpanStart`](@ref), [`SpanStop`](@ref). Empty gold span skipped; missing pred
vs present gold = 0.

[`default_metrics`](@ref) is allele + those four span metrics.

Panel ids: `{source}__assign={db_label}__sim={all|minus_holdout|holdout_only}`.
The assign FASTA is what tools receive (`assign=`). `sim_set` is which alleles
IgSim drew from the source FASTA. They can differ: e.g. sim `minus_holdout`
with assign = train FASTA (matched closed) vs assign = full FASTA (seen
reads, extra alleles in the gallery).

AIRR gold is the file's V/D/J calls after checking every token is a name in
the gold FASTA. A TSV annotated with a different database errors on load.

## Modes

| Mode | Use |
|------|-----|
| [`DiagnosticMode`](@ref) | Mid-training: small N, untimed by default, predictions off |
| [`FullReportMode`](@ref) | End report / CI: full N, one timed annotate pass, predictions + `summary.md` + `report.html` |
