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
        size_threshold_ignore=["apiref.md"], # a single autodocs page of every docstring is legitimately large
    ),
    plugins=[bib],
    pages=[
        "Introduction" => "index.md",
        "Getting Started" =>
            Any["General Usage" => "manual.md", "Scatterers and Boundaries" => "scatterers.md", "Quantities" => "quantities.md"],
        "Electromagnetics" => Any[
            "Scatterers" => Any["Spheres" => "electromagnetic/spheres.md"],
            "Excitations" => Any[
                "Plane Wave" => "electromagnetic/planeWave.md",
                "Dipoles" => "electromagnetic/dipoles.md",
                "Ring Currents" => "electromagnetic/ringCurrents.md",
                "Spherical Modes" => "electromagnetic/sphModes.md",
                "Uniform Static Field" => "electromagnetic/uniformStatic.md",
            ],
            "Radar Cross Section" => "electromagnetic/rcs.md",
        ],
        "Acoustics" => Any[
            "Scatterers" => Any["Sphere" => "acoustic/sphere.md", "Spheroid and Disc" => "acoustic/spheroid.md"],
            "Excitations" => Any["Plane Wave" => "acoustic/planeWave.md", "Monopole" => "acoustic/monopole.md"],
        ],
        "Numerical Details" => Any[
            "Coordinate Systems" => "numerics/coordinateSys.md",
            "Series for Spheres" => "numerics/sphereSeries.md",
            "Spheroidal Solution" => "numerics/spheroidSeries.md",
            "Units, Duality, and Rotations" => "numerics/details.md",
        ],
        "Examples" => Any["Code Verification" => "examples/verification.md", "Visualization of Fields" => "examples/visualization.md"],
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
