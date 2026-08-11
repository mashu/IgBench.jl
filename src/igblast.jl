# igblast.jl — IgBLAST.jl annotator (only built-in tool package dependency).

"""Annotate with NCBI IgBLAST via IgBLAST.jl."""
struct IgBLASTAnnotator{R} <: AbstractAnnotator
    name::String
    runner::R
    organism_param::String
end

function IgBLASTAnnotator(;
                          name::AbstractString = "igblast",
                          organism_param::AbstractString = "human",
                          additional_params::Dict{String,String} = Dict{String,String}(),
                          kwargs...)
    params = merge(Dict("organism" => String(organism_param),
                        "domain_system" => "imgt"),
                   additional_params)
    runner = IgBLAST.IgBLASTRunner(IgBLAST.IgBLASTn; additional_params = params, kwargs...)
    IgBLASTAnnotator(String(name), runner, String(organism_param))
end

tool_name(a::IgBLASTAnnotator) = a.name

function annotate(a::IgBLASTAnnotator, sequences, ids, germline::GermlinePaths; kwargs...)
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
        parsed = read_airr_calls(out_tsv; max_rows = nothing, skip_nonproductive = false)
        by_id = Dict(r.sequence_id => r for r in parsed)
        out = Vector{CallRecord}(undef, n)
        for i in 1:n
            r = get(by_id, String(ids[i]), nothing)
            if isnothing(r)
                out[i] = CallRecord(ids[i], sequences[i], "", "", "")
            else
                out[i] = r
            end
        end
        out
    end
end
