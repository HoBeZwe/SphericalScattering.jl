
@testset "Medium and Sphere" begin

    md = Medium(3.0, Float32(2.0))

    @test md isa Medium{Float64}

    mdf = Medium(Float32(2), Float32(1))

    @test mdf isa Medium{Float32}

    @test Medium{Float64}(mdf) isa Medium{Float64}

    pecsp = PECSphere(; radius=Float32(1.0))

    @test numlayers(pecsp) == 2

    ex_md = planeWave(; frequency=1e8, embedding=md)

    @test medium(pecsp, ex_md, 3.0) == ex_md.embedding

    @test medium(pecsp, ex_md, 0.5) == Medium(0.0, 0.0)

    @test SphericalScattering.wavenumber(pecsp, ex_md, 0.0) == 0.0

    sp = DielectricSphere(; radius=Float32(1.0), filling=md)

    @test medium(sp, ex_md, 0.5) == md
    @test medium(sp, ex_md, 1.5) == md
    @test medium(sp, planeWave(; frequency=1e8, embedding=mdf), 1.5) == mdf

    @test sp isa DielectricSphere{Float64,Float32}

    @test SphericalScattering.impedance(sp, planeWave(; frequency=1e8, embedding=mdf), 2.0) ≈ 0.707106 rtol = 1e-5

    @test SphericalScattering.wavenumber(sp, planeWave(; frequency=1e8, embedding=mdf), 2.0) ≈ 8.8857658e8 rtol = 1e-5

    md1 = Medium(15.0, -2.1) # innermost medium
    md2 = Medium(-10, -2.0) # innermost medium

    @test_throws ErrorException("Radii are not ordered ascendingly.") LayeredSphere(;
        radii=SVector(0.25, 0.5, 0.3), filling=SVector(md1, md2, md)
    )
    @test_throws ErrorException("Number of fillings does not match number of radii.") LayeredSphere(;
        radii=SVector(0.25, 0.3, 0.5), filling=SVector(md1, md2)
    )

    @test_throws ErrorException("Radii are not ordered ascendingly.") LayeredSpherePEC(;
        radii=SVector(0.25, 0.5, 0.3), filling=SVector(md1, md2)
    )
    @test_throws ErrorException("Number of fillings does not match number of radii.") LayeredSpherePEC(;
        radii=SVector(0.25, 0.3, 0.5), filling=SVector(md1)
    )

    spl = LayeredSphere(; radii=SVector(0.25, 0.3, 0.5), filling=SVector(md1, md2, md))

    @test numlayers(spl) == 4

    @test layer(spl, 0.1) == 1
    @test layer(spl, 0.24999) == 1
    @test layer(spl, 0.25) == 2
    @test layer(spl, 0.4) == 3
    @test layer(spl, 0.5) == 4
    @test layer(spl, 1.6) == 4

    ex = planeWave(; frequency=1e8, embedding=Medium(-3.0, 2.0))

    @test permittivity(spl, ex, 0.1) == 15.0
    @test permittivity(spl, ex, 0.24999) == 15.0
    @test permeability(spl, ex, 0.25) == -2.0
    @test permeability(spl, ex, 1.25) == 2.0

    spl = LayeredSpherePEC(; radii=SVector(0.25, 0.3, 0.5), filling=SVector(md1, md2))
    @test medium(spl, ex_md, 1.5) == md

    @test permittivity(spl, ex, 0.1) == 0.0
    @test permittivity(spl, ex, 0.25) == 15.0

    pts = [SVector(0.1, 0.0, 0.0), SVector(0.25, 0.0, 0.0)]

    @test permittivity(spl, ex, pts) == [0.0, 15.0]
    @test permeability(spl, ex, pts) == [0.0, -2.1]
end
