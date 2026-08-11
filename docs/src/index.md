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
# result.metrics is in-memory for logging; store optional (nothing / NullRunStore)
```

## Data

Describe cohorts with [`DatasetManifest`](@ref): [`SimSource`](@ref) (IgSim gold)
and [`AirrSource`](@ref) (real AIRR). `species` is a free string — never hardcoded.

Pass `cache = PanelCache("cache/run")` to [`run_suite`](@ref) so diagnostic steps
reuse the same frozen sequences.

## Modes

| Mode | Use |
|------|-----|
| [`DiagnosticMode`](@ref) | Mid-training: small N, light timing, predictions off by default |
| [`FullReportMode`](@ref) | End report / CI: full N, timings, predictions + `summary.md` |
