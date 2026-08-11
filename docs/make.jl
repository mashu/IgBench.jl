using Documenter
using IgBench

DocMeta.setdocmeta!(IgBench, :DocTestSetup, :(using IgBench); recursive = true)

makedocs(;
    modules = [IgBench],
    authors = "Mateusz Kaduk",
    sitename = "IgBench.jl",
    repo = Remotes.GitHub("mashu", "IgBench.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://mashu.github.io/IgBench.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Extending" => "extending.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/mashu/IgBench.jl",
    devbranch = "main",
    push_preview = true,
)
