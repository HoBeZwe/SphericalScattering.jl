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
        "Geometry" => Any["Coordinate System" => "coordinateSys.md", "Sphere Dimensions" => "scatterer.md"],
        "Excitations" => Any[
            "Plane Wave" => "planeWave.md",
            "Dipoles" => "dipoles.md",
            "Ring Currents" => "ringCurrents.md",
            "Spherical Modes" => "sphModes.md",
            "Uniform Static Field" => "uniformStatic.md",
        ],
        "Further Details" => "details.md",
        "Contributing" => "contributing.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/HoBeZwe/SphericalScattering.jl")
