# igblast.jl — IgBLAST annotator (implementation via Package Extension).

"""
Annotate with NCBI IgBLAST via IgBLAST.jl.

Requires IgBLAST.jl in the environment (weak dependency). Constructing without
IgBLAST loaded raises an error.
"""
struct IgBLASTAnnotator{R} <: AbstractAnnotator
    name::String
    runner::R
    organism_param::String
end

tool_name(a::IgBLASTAnnotator) = a.name

"""Build an [`IgBLASTAnnotator`](@ref); requires IgBLAST.jl to be loaded."""
function IgBLASTAnnotator(; kwargs...)
    igblast_annotator_impl(; kwargs...)
end

function igblast_annotator_impl(; kwargs...)
    error("IgBLAST.jl is not loaded. Add it to your environment " *
          "(Pkg.add(url=\"https://github.com/mashu/IgBLAST.jl\")) to use IgBLASTAnnotator.")
end

function annotate(::IgBLASTAnnotator, sequences, ids, germline::GermlinePaths; kwargs...)
    error("IgBLAST.jl is not loaded; cannot annotate with IgBLASTAnnotator")
end
