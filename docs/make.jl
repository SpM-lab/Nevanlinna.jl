using Nevanlinna
using Documenter

DocMeta.setdocmeta!(Nevanlinna, :DocTestSetup, :(using Nevanlinna); recursive=true)

makedocs(;
    modules=[Nevanlinna],
    authors="Hiroshi Shinaoka <h.shinaoka@gmail.com> and contributors",
    repo="https://github.com/shinaoka/Nevanlinna.jl/blob/{commit}{path}#{line}",
    sitename="Nevanlinna.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://shinaoka.github.io/Nevanlinna.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/shinaoka/Nevanlinna.jl",
    devbranch="main",
)
