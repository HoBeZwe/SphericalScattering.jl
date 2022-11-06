using SphericalScattering
using Test

using JuliaFormatter
using StaticArrays
using BEAST
using CompScienceMeshes
using LinearAlgebra


ϑ = range(0.0 + 0.001; stop=π - 0.0001, length=18)  # 10° steps
ϕ = range(0.0; stop=2π, length=36) # 10° steps

function getDefaultPoints(r::Float64)
    points_cart = [point(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]       # convert points to cartesian components
    points_sph  = [point(r, θ, φ) for θ in ϑ, φ in ϕ]

    return points_cart, points_sph
end


@testset "Testing SphericalScattering functionality" begin
    @testset "Test dipoles" begin
        include("dipoles.jl")
    end
    @testset "Test plane waves" begin
        include("planeWave.jl")
    end
    @testset "Test ring currents" begin
        include("ringCurrents.jl")
    end
    @testset "Test spherical modes" begin
        include("sphericalModes.jl")
    end
    @testset "Test uniform field" begin
        include("uniformField.jl")
    end
    @testset "Test formatting of files" begin
        pkgpath = pkgdir(SphericalScattering)   # path of this package including name
        @test format(pkgpath)                   # check whether files are formatted according to the .JuliaFormatter.toml 
    end
end
