using Documenter
using FluxionFields


makedocs(
    modules = [FluxionFields],
    format = :html,
    sitename = "FluxionFields.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/fluxion/FluxionFields.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
