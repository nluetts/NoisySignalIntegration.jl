using Documenter
using MCIntegrate
using Pkg

Pkg.activate(".") # activate the docs environment

# include doctests from docstrings
DocMeta.setdocmeta!(MCIntegrate, :DocTestSetup, :(using MCIntegrate); recursive=true)

makedocs(
    sitename="NoisySignalIntegration.jl Documentation",
    modules=[MCIntegrate],
    pages = [
        "Home" => "index.md",
        "Usage Guide" => "guide.md",
        "API reference" => "api_ref.md"
    ],
    highlightsig=true
)