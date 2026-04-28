using Documenter
using DiscreteCurves

DocMeta.setdocmeta!(DiscreteCurves, :DocTestSetup,
    :(using DiscreteCurves, StaticArrays); recursive=true)

makedocs(;
    modules  = [DiscreteCurves],
    sitename = "DiscreteCurves.jl",
    authors  = "Cornelia Vizman, Sava Girban",
    format   = Documenter.HTML(;
        prettyurls    = get(ENV, "CI", nothing) == "true",
        canonical     = "https://sgirban.github.io/DiscreteCurves.jl",
        edit_link     = "master",
        assets        = String[],
    ),
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Mathematics"     => "mathematics.md",
    ],
    checkdocs = :exports,
    warnonly  = true,
)
