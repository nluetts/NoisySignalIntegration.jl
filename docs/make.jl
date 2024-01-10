import Pkg
Pkg.activate(@__DIR__)

push!(LOAD_PATH, "../src")

using Documenter
using Plots # this needs to be imported here so that Require.jl includes animate_draws when building docs
using NoisySignalIntegration

# include doctests from docstrings
DocMeta.setdocmeta!(NoisySignalIntegration, :DocTestSetup, :(using NoisySignalIntegration); recursive=true)

makedocs(
    sitename="NoisySignalIntegration.jl",
    modules=[NoisySignalIntegration],
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Package Overview" => "overview.md",
        "Usage Guide" => "guide.md",
        "Examples" => "examples.md",
        "Baseline Handling" => "baseline.md",
        "Internals" => "internals.md",
        "API reference" => "API.md"
    ],
    highlightsig=true,
    checkdocs=:none
)
