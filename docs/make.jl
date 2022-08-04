using QRDecoders
using Documenter

DocMeta.setdocmeta!(QRDecoders, :DocTestSetup, :(using QRDecoders); recursive=true)

makedocs(;
    modules=[QRDecoders],
    authors="rex <1073853456@qq.com> and contributors",
    repo="https://github.com/RexWzh/QRDecoders.jl/blob/{commit}{path}#{line}",
    sitename="QRDecoders.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RexWzh.github.io/QRDecoders.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RexWzh/QRDecoders.jl",
    devbranch="main",
)
