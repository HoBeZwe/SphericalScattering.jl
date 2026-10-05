
@testitem "Medium and Sphere" setup = [Setup] begin

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


@testitem "Scatterer hierarchy" setup = [Setup] begin

    SS = SphericalScattering

    @testset "Geometry and condition" begin

        # --- the three ways of constructing a sphere with a parameter-free condition coincide
        sp = HardSphere(; radius=1.0)

        @test sp === SS.Sphere{SoundHard}(; radius=1.0)
        @test sp === SS.Sphere(; radius=1.0, boundary=SoundHard())
        @test SoftSphere(; radius=1.0) === SS.Sphere{SoundSoft}(; radius=1.0)

        @test sp isa HardSphere
        @test sp isa Scatterer{SoundHard}
        @test sp isa Scatterer{<:AcousticBoundary}
        @test !(sp isa Scatterer{<:ElectromagneticBoundary})

        # --- the condition is a type parameter and a field, the former costing no memory for a tag
        @test sp.boundary === SoundHard()
        @test sizeof(sp) == sizeof(1.0)

        # --- likewise for the electromagnetic spheres, whose positional constructors are kept
        pec = PECSphere(; radius=1.0)

        @test pec === SS.Sphere{PEC}(; radius=1.0)
        @test pec === SS.Sphere(; radius=1.0, boundary=PEC())
        @test pec === PECSphere(1.0)
        @test pec isa Scatterer{PEC}
        @test sizeof(pec) == sizeof(1.0)

        # --- a condition requiring data is passed as a value; the alias keeps the parameters of the former type
        md = Medium(2.0, 1.0)
        dsp = DielectricSphere(; radius=Float32(1.0), filling=md)

        @test dsp == SS.Sphere(; radius=Float32(1.0), boundary=Dielectric(md))
        @test dsp == DielectricSphere(Float32(1.0), md)
        @test dsp isa SS.Sphere{<:Dielectric}
        @test dsp isa DielectricSphere{Float64,Float32}
        @test dsp.boundary.filling == md

        # --- a condition which is not parameter-free cannot be given as a type
        @test_throws ErrorException SS.Sphere{Dielectric}(; radius=1.0)
        @test_throws ErrorException SS.Sphere{ElectromagneticBoundary}(; radius=1.0)

        # --- a layered sphere stores its outermost radius as the radius of the sphere, and the inner interfaces,
        #     the fillings of the shells and the core as its condition
        md1, md2 = Medium(15.0, -2.1), Medium(-10.0, -2.0)

        spl = LayeredSphere(; radii=SVector(0.25, 0.5, 1.0), filling=SVector(md1, md2, md))

        @test spl == LayeredSphere(SVector(0.25, 0.5, 1.0), SVector(md1, md2, md))
        @test spl isa SS.Sphere{<:Layered{<:Dielectric}}
        @test !(spl isa LayeredSpherePEC)
        @test spl.radius == 1.0
        @test spl.boundary.radii == SVector(0.25, 0.5)
        @test spl.boundary.filling == SVector(md2, md)
        @test spl.boundary.core == Dielectric(md1)
        @test SS.layerRadii(spl) == SVector(0.25, 0.5, 1.0)  # the views the algorithms work with
        @test SS.layerFillings(spl) == SVector(md1, md2, md)

        spp = LayeredSpherePEC(; radii=SVector(0.25, 0.5, 1.0), filling=SVector(md1, md2))

        @test spp == LayeredSpherePEC(SVector(0.25, 0.5, 1.0), SVector(md1, md2))
        @test spp isa SS.Sphere{<:Layered{PEC}}
        @test !(spp isa LayeredSphere)
        @test spp.boundary.core === PEC()
        @test SS.layerRadii(spp) == SVector(0.25, 0.5, 1.0)
        @test SS.layerFillings(spp) == SVector(md1, md2)

        # --- without shells, the core remains
        @test LayeredSpherePEC(; radii=SVector(1.0), filling=SVector{0,Medium{Float64}}()).boundary.radii == SVector{0,Float64}()

        # --- the thin impedance layer is a condition as well; the alias keeps the parameters of the former type
        spj = DielectricSphereThinImpedanceLayer(; radius=1.0, thickness=0.01, thinlayer=md2, filling=md)

        @test spj isa SS.Sphere{<:ThinImpedanceLayer}
        @test spj isa DielectricSphereThinImpedanceLayer{Float64,Float64}
        @test spj.boundary.thickness == 0.01
        @test spj.boundary.thinlayer == md2
        @test DielectricSphereThinImpedanceLayer(1, Float32(0.01), md2, md) isa DielectricSphereThinImpedanceLayer{Float32,Float64}

        # --- the interior is a question of geometry alone
        for sp in (HardSphere(; radius=1.0), SoftSphere(; radius=1.0), pec, DielectricSphere(; radius=1.0, filling=md), spl, spp, spj)
            @test SS.isinside(sp, SVector(0.5, 0.0, 0.0))
            @test !SS.isinside(sp, SVector(1.5, 0.0, 0.0))
        end

        # --- a spheroid is a scatterer of its own, not a sphere
        sph = Spheroid{SoundSoft}(; equatorialRadius=2.0, polarRadius=1.0)

        @test sph isa Scatterer{SoundSoft}
        @test !(sph isa SS.Sphere)
        @test Disc(SoundHard; radius=1.0) isa Scatterer{SoundHard}

        # --- the shapes are subtypes of `Spheroid`, the condition being a field as for a sphere
        @test sph isa OblateSpheroid{SoundSoft}
        @test sph isa Spheroid{SoundSoft}
        @test sph.boundary === SoundSoft()
        @test sizeof(sph) == 5 * sizeof(1.0) # the semifocal distance, ξ₀ and the axis: the tag costs nothing
        @test Disc(SoundHard; radius=1.0) isa OblateSpheroid{SoundHard}
    end

    @testset "Physics" begin

        md = Medium(2.0, 1.0)

        for sp in (
            PECSphere(; radius=1.0),
            DielectricSphere(; radius=1.0, filling=md),
            LayeredSphere(; radii=SVector(0.5, 1.0), filling=SVector(md, md)),
            LayeredSpherePEC(; radii=SVector(0.5, 1.0), filling=SVector(md)),
            DielectricSphereThinImpedanceLayer(; radius=1.0, thickness=0.01, thinlayer=md, filling=md),
        )
            @test sp isa SS.Sphere # every scatterer of a spherical geometry is a `Sphere` now
            @test sp isa Scatterer{<:ElectromagneticBoundary}
            @test !(sp isa Scatterer{<:AcousticBoundary})
        end

        for ex in
            (planeWave(; frequency=f), HertzianDipole(; frequency=f, position=SVector(0.0, 0.0, 3.0)), SphericalModeTE(; frequency=f))
            @test ex isa SS.ElectromagneticExcitation
        end

        for ex in (SS.Acoustic.planeWave(; frequency=f), SS.Acoustic.monopole(; frequency=f, position=SVector(0.0, 0.0, 3.0)))
            @test ex isa SS.AcousticExcitation
            @test !(ex isa SS.ElectromagneticExcitation)
        end
    end

    @testset "Scatterer and excitation of different physics" begin

        # the mismatch is rejected by dispatch, before any parallel loop: the error must not depend on
        # the number of threads, as it would if thrown inside a loop
        errEM = ErrorException(
            "An electromagnetic excitation requires a scatterer with an electromagnetic boundary condition, such as a `PECSphere`."
        )
        errAC = ErrorException(
            "An acoustic excitation requires a scatterer with an acoustic boundary condition, such as a `HardSphere`, a `SoftSphere` or a `Spheroid`.",
        )

        points = [SVector(2.0, 0.0, 0.0), SVector(0.0, 2.0, 0.0)]

        for sp in (HardSphere(; radius=1.0), Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=0.5))
            for ex in (planeWave(; frequency=f), HertzianDipole(; frequency=f, position=SVector(0.0, 0.0, 3.0)))
                @test_throws errEM scatteredfield(sp, ex, ElectricField(points))
                @test_throws errEM field(sp, ex, ElectricField(points))
            end
        end

        pec = PECSphere(; radius=1.0)

        for ex in (SS.Acoustic.planeWave(; frequency=f), SS.Acoustic.monopole(; frequency=f, position=SVector(0.0, 0.0, 3.0)))
            for quantity in (Pressure, PressureTrace, PressureNormalGradient, PressureJump)
                @test_throws errAC scatteredfield(pec, ex, quantity(points))
            end
            @test_throws errAC field(pec, ex, Pressure(points))
        end
    end
end
