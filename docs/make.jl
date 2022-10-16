using Documenter
using QRDecoders

DocMeta.setdocmeta!(QRDecoders, :DocTestSetup, :(using QRDecoders); recursive=true)

makedocs(;
    modules=[QRDecoders],
    sitename="QRDecoders.jl"
)

deploydocs(;
    repo="github.com/JuliaImages/QRDecoders.jl",
)
