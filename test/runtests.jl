using SphericalScattering
using Test

using JuliaFormatter
using StaticArrays

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
