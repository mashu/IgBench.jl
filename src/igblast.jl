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

# Extension hook — declare only (no methods). IgBenchIgBLASTExt adds the method.
# Defining a stub method here and overwriting it in the ext breaks precompilation
# on Julia 1.12+ ("Method overwriting is not permitted").
function igblast_annotator_impl end

"""Build an [`IgBLASTAnnotator`](@ref); requires IgBLAST.jl to be loaded."""
function IgBLASTAnnotator(; kwargs...)
    isempty(methods(igblast_annotator_impl)) &&
        error("IgBLAST.jl is not loaded. Add it to your environment " *
              "(Pkg.add(url=\"https://github.com/mashu/IgBLAST.jl\")) to use IgBLASTAnnotator.")
    igblast_annotator_impl(; kwargs...)
end

# `annotate(::IgBLASTAnnotator, ...)` is provided by IgBenchIgBLASTExt only.
