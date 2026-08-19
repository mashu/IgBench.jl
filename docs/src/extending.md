# Extending

## New tool (IgFormer example)

IgBench does **not** depend on IgFormer. In the IgFormer project:

```julia
using IgBench

igf = CallableAnnotator("igformer") do seqs, ids, germline
    # annotate against `germline` (assign FASTA paths); return Vector{CallRecord}
    CallRecord[]  # replace
end

suite = suite_from_manifest(manifest, [IgBLASTAnnotator(; organism_param="human"), igf])
run_suite(suite; mode = FullReportMode(), store = DirectoryRunStore("runs/full/end"))
```

SwiftIG (no Julia package):

```julia
SwiftIGAnnotator(; bin = get(ENV, "SWIFTIG_BIN", "swiftig"), threads = 8)
```

## New metric

```julia
struct MyMetric <: AbstractMetric end
IgBench.metric_name(::MyMetric) = "mine"
function IgBench.evaluate(::MyMetric, pred, ref)
    MetricValue(v=0.0, d=0.0, j=0.0, n=length(pred))
end
```

## New panel

Subtype [`AbstractPanel`](@ref) and implement `load_panel` → [`PanelData`](@ref).
