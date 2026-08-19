module IgBenchIgBLASTExt

using IgBench
using IgBLAST

# Add a method to the parent hook (parent declares the function with no methods).
function IgBench.igblast_annotator_impl(;
                                        name::AbstractString = "igblast",
                                        organism_param::AbstractString = "human",
                                        num_threads::Integer = 8,
                                        additional_params::Dict{String,String} = Dict{String,String}(),
                                        kwargs...)
    params = merge(Dict("organism" => String(organism_param),
                        "domain_system" => "imgt"),
                   additional_params)
    runner = IgBLAST.IgBLASTRunner(IgBLAST.IgBLASTn;
                                   additional_params = params,
                                   num_threads,
                                   kwargs...)
    IgBench.IgBLASTAnnotator(String(name), runner, String(organism_param))
end

function IgBench.annotate(a::IgBench.IgBLASTAnnotator, sequences, ids,
                          germline::IgBench.GermlinePaths; kwargs...)
    n = length(ids)
    length(sequences) == n || error("sequences/ids length mismatch")
    isnothing(germline.d) && error("IgBLASTAnnotator requires a D germline FASTA path")
    mktempdir() do dir
        query = joinpath(dir, "query.fasta")
        out_tsv = joinpath(dir, "out.tsv")
        open(query, "w") do io
            for (id, seq) in zip(ids, sequences)
                println(io, '>', id)
                println(io, seq)
            end
        end
        dbs = IgBLAST.VDJGermlines(germline.v, germline.d, germline.j)
        a.runner(query, dbs, out_tsv)
        parsed = IgBench.read_airr_calls(out_tsv; max_rows = nothing, skip_nonproductive = false)
        by_id = Dict(r.sequence_id => r for r in parsed)
        out = Vector{IgBench.CallRecord}(undef, n)
        for i in 1:n
            r = get(by_id, String(ids[i]), nothing)
            if isnothing(r)
                out[i] = IgBench.CallRecord(ids[i], sequences[i], "", "", "")
            else
                out[i] = r
            end
        end
        out
    end
end

end # module
