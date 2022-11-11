using SphericalScattering
using Documenter

DocMeta.setdocmeta!(SphericalScattering, :DocTestSetup, :(using SphericalScattering); recursive=true)

makedocs(;
    modules=[SphericalScattering],
    authors="Bernd Hofmann <Bernd.Hofmann@tum.de> and contributors",
    repo="https://github.com/HoBeZwe/SphericalScattering.jl/tree/master",
    sitename="SphericalScattering.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://HoBeZwe.github.io/SphericalScattering.jl", assets=String[]
    ),
    pages=[
        "Introduction" => "index.md",
        "Manual" => "manual.md",
        "Implementation Details" => Any[
            "General" => "general.md",
            "Plane Wave" => "planeWave.md",
            "Ring Currents" => "ringCurrents.md",
            "Dipoles" => "dipoles.md",
            "Spherical Modes" => "sphModes.md",
        ],
        "Contributing" => "contributing.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/HoBeZwe/SphericalScattering.jl")
