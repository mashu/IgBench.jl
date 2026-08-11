# IgBench.jl

[![Build Status](https://github.com/mashu/IgBench.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBench.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBench.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBench.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/IgBench.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBench.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Extensible V(D)J benchmark harness: pluggable annotators, metrics, and
species-agnostic panels. Library-first (`run_suite`) for embedding from
IgFormer; optional standalone script for full reports.

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
src = SimSource(; id="sim", germline=gp, species="human", n=256, seed=1,
                holdout_v=["IGHV1-69*01"])
man = DatasetManifest("demo"; sim=[src], airr=AirrSource[])

igblast = IgBLASTAnnotator(; organism_param="human")
# IgFormer (or any model) — no IgBench → IgFormer dependency:
model = CallableAnnotator("igformer") do seqs, ids, germline
    # ... annotate_many → Vector{CallRecord}
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

```
runs/<mode>/<id>/
  manifest.json
  panels.json
  predictions/<panel>__<tool>.airr.tsv.gz
  metrics.jsonl          # one row per panel×compare×metric (+ step tags)
  metrics.json
  timing.json
  summary.md             # full mode
```

## Extending

- **Tool:** subtype `AbstractAnnotator` or use `CallableAnnotator` / `SwiftIGAnnotator`
- **Metric:** subtype `AbstractMetric`, implement `evaluate`
- **Panel:** subtype `AbstractPanel`, implement `load_panel`

## License

MIT — see [LICENSE](LICENSE).
