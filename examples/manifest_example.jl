# Example DatasetManifest — edit paths for your machine.
# Species tags are free strings (not hardcoded to human/rhesus).
#
#   julia --project=. scripts/run_bench.jl --manifest examples/manifest_example.jl --mode diagnostic --out /tmp/igbench

using IgBench

# Replace with your germline / AIRR paths:
gp = GermlinePaths(;
    v = joinpath(@__DIR__, "..", "test", "fixtures", "V.fasta"),
    d = joinpath(@__DIR__, "..", "test", "fixtures", "D.fasta"),
    j = joinpath(@__DIR__, "..", "test", "fixtures", "J.fasta"),
)

DatasetManifest("example";
    sim = [
        SimSource(; id = "sim_toy", germline = gp, species = "toy",
                  n = 32, seed = 1, holdout_v = String[]),
    ],
    airr = AirrSource[],
)
