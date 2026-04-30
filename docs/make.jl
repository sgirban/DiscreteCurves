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
        sidebar_sitename = true,
    ),
    pages = [
        "Home"          => "index.md",
        "Getting Started" => "getting_started.md",
        "Mathematics"   => "mathematics.md",
        "API Reference" => [
            "Constructors"          => "api/constructors.md",
            "Topology & Access"     => "api/topology.md",
            "Lengths"               => "api/lengths.md",
            "Curvature"             => "api/curvature.md",
            "Tangents & Normals"    => "api/vectors.md",
            "Geometric Properties"  => "api/geometric_properties.md",
            "Orientation"           => "api/orientation.md",
            "Data Attachments"      => "api/data.md",
            "Curve Generators"      => "api/generators.md",
            "Flows"                 => "api/flows.md",
            # "Macros"                => "api/macros.md",
            # "Plotting"              => "api/plotting.md",
        ],
    ],
    checkdocs = :exports,
    warnonly  = true,
)
