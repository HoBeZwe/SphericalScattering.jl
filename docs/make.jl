using SphericalScattering
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(;
    modules=[SphericalScattering],
    authors="Bernd Hofmann <Bernd.Hofmann@tum.de> and contributors",
    sitename="SphericalScattering.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HoBeZwe.github.io/SphericalScattering.jl",
        edit_link="master",
        assets=String[],
        collapselevel=2,
        sidebar_sitename=true,
    ),
    plugins=[bib],
    pages=[
        "Introduction" => "index.md",
        "Manual" => Any["General Usage" => "manual.md", "Application Examples" => "application.md"],
        "Geometry" => Any["Coordinate System" => "coordinateSys.md", "Sphere Dimensions" => "scatterer.md"],
        "Excitations" => Any[
            "Electromagnetic" => Any[
                "Plane Wave" => "electromagnetic/planeWave.md",
                "Dipoles" => "electromagnetic/dipoles.md",
                "Ring Currents" => "electromagnetic/ringCurrents.md",
                "Spherical Modes" => "electromagnetic/sphModes.md",
                "Uniform Static Field" => "electromagnetic/uniformStatic.md",
            ],
            "Acoustic" => Any["Plane Wave" => "acoustic/planeWave.md", "Monopole" => "acoustic/monopole.md"],
        ],
        "Further Details" => "details.md",
        "Contributing" => "contributing.md",
        "References" => "references.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(;
    repo="github.com/HoBeZwe/SphericalScattering.jl",
    target="build",
    push_preview=true,
    forcepush=true,
    versions=["stable" => "v^", "v#.#", "v0.5.0", "v0.4.0", "v0.3.0", "v0.2.0", "v0.1.2", "v0.1.1", "dev" => "dev"],
)
