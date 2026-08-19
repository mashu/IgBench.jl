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

const IGBLAST_PKG = Base.PkgId(
    Base.UUID("0b79348f-fd03-46c5-b933-bdb7f1bed3f3"), "IgBLAST")

igblast_extension_ready() = !isempty(methods(igblast_annotator_impl))

"""Load IgBLAST.jl if it is installed so the package extension can activate."""
function load_igblast_extension!()
    igblast_extension_ready() && return
    isnothing(Base.locate_package(IGBLAST_PKG)) && return
    Base.require(IGBLAST_PKG)
end

"""Build an [`IgBLASTAnnotator`](@ref); requires IgBLAST.jl in the environment."""
function IgBLASTAnnotator(; kwargs...)
    load_igblast_extension!()
    igblast_extension_ready() ||
        error("IgBLAST.jl is not loaded. Add it to your environment " *
              "(Pkg.add(url=\"https://github.com/mashu/IgBLAST.jl\")) to use IgBLASTAnnotator.")
    Base.invokelatest(igblast_annotator_impl; kwargs...)
end

# `annotate(::IgBLASTAnnotator, ...)` is provided by IgBenchIgBLASTExt only.
