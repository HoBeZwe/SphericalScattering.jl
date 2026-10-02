using TestItemRunner

using SphericalScattering
using Test

using JuliaFormatter
using StaticArrays
using BEAST
using CompScienceMeshes
using LinearAlgebra


# ----- points on spherical grid
ϑ = range(0.0; stop=π, length=18)  # 10° steps
ϕ = range(0.0; stop=2π, length=36) # 10° steps

function getDefaultPoints(r::Float64)
    points_cart = [point(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]       # convert points to cartesian components
    points_sph = [point(r, θ, φ) for θ in ϑ, φ in ϕ]

    return points_cart, points_sph
end


# ----- interface to BEAST
function (lc::Excitation)(p)

    F_cart = field(lc, ElectricField([p]))

    return F_cart[1]
end

BEAST.cross(::BEAST.NormalVector, p::Excitation) = CrossTraceMW(p)
BEAST.scalartype(p::Excitation) = Complex{typeof(p.embedding.ε)}


# ----- variables used in all tests
spRadius = 1.0 # radius of sphere

Γ = meshsphere(spRadius, 0.45)
RT = raviartthomas(Γ)

𝜇 = SphericalScattering.μ0
𝜀 = SphericalScattering.ε0

c = 1 / sqrt(𝜇 * 𝜀)

points_cartFF, points_sphFF = getDefaultPoints(1.0)
points_cartNF, points_sphNF = getDefaultPoints(5.0)
points_cartNF_inside, ~ = getDefaultPoints(0.5)

# ----- testsets
@testset verbose = true "Testing SphericalScattering functionality" begin

    @testset "Medium and Sphere" begin
        include("sphere.jl")
    end

    @testset "Test dipoles" begin
        include("electromagnetics/dipoles.jl")
    end

    @testset "Test plane waves" begin
        include("electromagnetics/planeWave.jl")
        include("electromagnetics/planeWave_dielectric.jl")
    end

    @testset "Test ring currents" begin
        include("electromagnetics/ringCurrents.jl")
    end

    @testset "Test spherical modes" begin
        include("electromagnetics/sphericalModes.jl")
    end

    @testset "Test uniform field" begin
        include("electromagnetics/uniformField.jl")
    end

    @testset "Test acoustics" begin
        include("acoustics/planeWave.jl")
        include("acoustics/monopole.jl")
        include("acoustics/spheroid.jl")
    end

    @testset "Test coordinate transforms" begin
        include("coordinateTransforms.jl")
    end

    @testset "Test utils" begin
        include("utils.jl")
    end

    @testset "Test formatting of files" begin
        pkgpath = pkgdir(SphericalScattering)   # path of this package including name
        @test format(pkgpath, overwrite=false)  # check whether files are formatted according to the .JuliaFormatter.toml 
    end
end
