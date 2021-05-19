import Pkg
Pkg.activate(@__DIR__)

using Documenter
using NoisySignalIntegration

# include doctests from docstrings
DocMeta.setdocmeta!(NoisySignalIntegration, :DocTestSetup, :(using NoisySignalIntegration); recursive=true)

makedocs(
    sitename="NoisySignalIntegration.jl",
    modules=[NoisySignalIntegration],
    pages = [
        "Home" => "index.md",
        "Package Overview" => "overview.md",
        "Usage Guide" => "guide.md",
        "Examples" => "examples.md",
        "Internals" => "internals.md",
        "API reference" => "API.md"
    ],
    highlightsig=true
)