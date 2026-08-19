# IgBench.jl

[![Build Status](https://github.com/mashu/IgBench.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBench.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBench.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBench.jl/badge.svg)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/IgBench.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBench.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Extensible V(D)J benchmark harness: pluggable annotators, metrics, and
species-agnostic panels. Library-first (`run_suite`) for embedding from
IgFormer; optional standalone script for full reports.

### Published metrics

**Allele calling** is [`AlleleAccuracy`](@ref) (name `"allele"`):

1. Empty / `NA` / `.` **gold** → skip (not in the denominator) on V, D, and J.
2. Present gold + empty **pred** → score `0` (false negative).
3. Full-string match → `1`. Else split on **commas only** (slash is part of an
   IMGT dual name such as `IGHV3-23*01/IGHV3-23D*01`); any pred token equal to
   any gold token → **`1 / n_pred`**.

**Segmentation** is independent of the allele score:

| name | meaning |
|---|---|
| `span_iou` | intersection over union |
| `span_exact` | start **and** stop equal gold |
| `span_start` | start equals gold |
| `span_stop` | stop equals gold |

Empty gold span is skipped; missing pred vs present gold scores `0`.

**Timing** is wall-clock seconds of the **scored** annotate pass, plus `n` and
reads/s. Tools are not annotated twice.

Panel ids name the **assign FASTA** and, for IgSim, the **sim allele-set**:
`{source}__assign={db}__sim={all|minus_holdout|holdout_only}`.

**Built-in tool:** [IgBLAST.jl](https://github.com/mashu/IgBLAST.jl) via a Package
Extension (optional at load time). SwiftIG (CLI) and IgFormer plug in with **no**
package dependency.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgBench.jl")
```

## Quick start (library)

```julia
using IgBench, IgSim

gp = GermlinePaths(; v="V.fasta", d="D.fasta", j="J.fasta")
src = SimSource(; id="sim", db_label="KI+1KGP", germline=gp, species="human",
                n=256, seed=1, holdout_v=["IGHV1-69*01"])
man = DatasetManifest("demo"; sim=[src], airr=AirrSource[])

igblast = IgBLASTAnnotator(; organism_param="human")
model = CallableAnnotator("igformer") do seqs, ids, germline
    # annotate against germline (assign FASTA); return Vector{CallRecord}
    CallRecord[]
end

suite = suite_from_manifest(man, [igblast, model])
result = run_suite(suite;
                   mode = DiagnosticMode(max_sequences=128),
                   store = DirectoryRunStore("runs/diagnostic/1000"),
                   step = 1000, tags = (; epoch = 3))
```

## Standalone

```bash
julia --project=. scripts/run_bench.jl \
  --manifest examples/manifest_example.jl \
  --mode full \
  --out runs/full/demo
```

## Run artifacts (dashboard-ready)

A run is a **directory**. Filenames inside it are fixed (`report.html`, `summary.md`,
`metrics.json`, …). A new report is a new `--out` directory; the same path overwrites.

```
runs/<mode>/<id>/
  manifest.json
  panels.json
  predictions/<panel>__<tool>.airr.tsv.gz
  metrics.jsonl          # one row per panel×compare×metric (+ step tags)
  metrics.json
  timing.json
  summary.md             # full mode
  report.html            # full mode; self-contained, open in a browser
```

## Extending

- **Tool:** subtype `AbstractAnnotator` or use `CallableAnnotator` / `SwiftIGAnnotator`
- **Metric:** subtype `AbstractMetric`, implement `evaluate`
- **Panel:** subtype `AbstractPanel`, implement `load_panel`

## License

MIT — see [LICENSE](LICENSE).
