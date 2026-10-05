
# The prolate spheroid, with the checks of `oblate.jl` for its shape; the coordinates are tested in
# `geometry/coordinateTransforms.jl`, the machinery independent of the shape in `machinery.jl`.

@testitem "Prolate spheroid" setup = [Setup] begin

    SS = SphericalScattering

    @testset "Types" begin

        # --- the shape follows from the radii
        sp = Spheroid{SoundHard}(; equatorialRadius=0.6, polarRadius=1.0)

        @test sp isa ProlateSpheroid{SoundHard,Float64}
        @test sp isa Spheroid{SoundHard}
        @test Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=0.6) isa OblateSpheroid{SoundHard,Float64}

        @test sp.semifocal ≈ 0.8 rtol = 1e-15
        @test equatorialRadius(sp) ≈ 0.6 rtol = 1e-14
        @test polarRadius(sp) ≈ 1.0 rtol = 1e-14
        @test !isdisc(sp)
        @test sp.boundary === SoundHard()

        # --- the sphere, the needle and an oblate shape are rejected
        @test_throws ErrorException Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=1.0)
        @test_throws ErrorException ProlateSpheroid{SoundHard}(; equatorialRadius=0.0, polarRadius=1.0)
        @test_throws ErrorException ProlateSpheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=0.6)

        # --- the interior
        @test SS.isinside(sp, SVector(0.0, 0.0, 0.9))
        @test SS.isinside(sp, SVector(0.5, 0.0, 0.0))
        @test !SS.isinside(sp, SVector(0.0, 0.0, 1.1))
        @test !SS.isinside(sp, SVector(0.65, 0.0, 0.0))
    end

    @testset "Wave functions" begin

        κ = 1.7

        # --- the Wronskian of the two real radial solutions is 1 / (c (ξ² - 1)), which fixes the radial coordinate
        for (m, n, ξ) in ((0, 0, 1.3), (1, 3, 1.05), (2, 5, 2.4))
            R₁ = SS.rmn(m, n, κ, [ξ]; spheroid=:prolate, kind=1)
            R₂ = SS.rmn(m, n, κ, [ξ]; spheroid=:prolate, kind=2)

            W = only(R₁.value) * only(R₂.derivative) - only(R₂.value) * only(R₁.derivative)

            @test real(W) ≈ 1 / (κ * (ξ^2 - 1)) rtol = 1e-10
        end

        # --- the angular functions reduce to the Legendre functions as c → 0
        for (m, n) in ((0, 2), (1, 3), (2, 4))
            @test SS.spheroidalAngular(:prolate, m, n, 1e-7, 0.37).value ≈ SS.Plm(0.37, n, m) rtol = 1e-6
        end

        # --- the outgoing kind behaves as jⁿ⁺¹ exp(-j c ξ) / (c ξ), as for the oblate shape; kind 3 is incoming
        ξs = (20.0, 30.0, 40.0)
        Q(kind) = [only(SS.rmn(0, 1, κ, [ξ]; spheroid=:prolate, kind=kind).value) * κ * ξ * cis(κ * ξ) for ξ in ξs]

        @test maximum(abs.(Q(4) .- im^2)) < 2e-2
        @test maximum(abs.(Q(3) .- Q(3)[1])) > 1.0
    end
end


@testitem "Prolate scattering" setup = [Setup] begin

    SS = SphericalScattering

    freq(k) = k * c / (2π)
    quiet(f) = redirect_stdout(f, devnull) # the check of the truncation may report

    locate(sp, η, φ) = SS.frame(sp) * SS.cartesianCoordinates(sp, SVector(sp.ξ₀, η, φ))
    grid = [(η, φ) for η in (-0.9, -0.4, 0.1, 0.6, 0.95) for φ in (0.3, 2.2, 4.4)]

    planeWave = SS.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))
    monopole = SS.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.5, -1.0, 2.0))

    @testset "Boundary conditions" begin

        # tilted, a moderate and a thin spheroid, the latter close to the excluded line segment
        for BC in (SoundSoft, SoundHard), (a, axis) in ((0.6, SVector(0.2, -0.3, 1.0)), (0.05, SVector(1.0, 0.4, 0.2)))
            sp = Spheroid{BC}(; equatorialRadius=a, polarRadius=1.0, axis=axis)
            points = [locate(sp, η, φ) for (η, φ) in grid]

            for ex in (planeWave, monopole)
                md = quiet(() -> SS.modes(sp, ex))

                q = BC === SoundSoft ? PressureTrace(points) : PressureNormalGradient(sp, points)

                @test maximum(abs.(field(sp, ex, md, q))) / maximum(abs.(field(ex, q))) < (BC === SoundSoft ? 1e-11 : 1e-9)
            end
        end
    end

    @testset "Sphere bridge" begin

        # As the polar radius approaches the equatorial one, the prolate spheroid degenerates into a sphere, the
        # deviation from the independently validated series of the sphere being of the order of the eccentricity
        ex = SS.Acoustic.planeWave(; frequency=freq(2.0), direction=SVector(0.0, 0.0, 1.0))
        points = [SVector(0.0, 0.0, 3.0), SVector(2.0, 1.0, -1.5), SVector(-3.0, 0.0, 0.0)]

        for (BC, SphereType) in ((SoundSoft, SoftSphere), (SoundHard, HardSphere))
            reference = [scatteredfield(SphereType(; radius=1.0), ex, Pressure([p]))[1] for p in points]

            deviation = map((1e-2, 1e-4)) do δ
                sp = ProlateSpheroid{BC}(; equatorialRadius=1.0, polarRadius=1.0 + δ)
                md = quiet(() -> SS.modes(sp, ex; M=2, N=16))

                maximum(abs.([scatteredfield(sp, ex, md, p, Pressure(nothing)) for p in points] .- reference) ./ abs.(reference))
            end

            @test deviation[2] < deviation[1] / 50
            @test deviation[2] < 1e-3
        end
    end

    @testset "Radiation condition" begin

        # p r exp(j k r) approaches a constant for the outgoing radial function, which no boundary condition would
        # distinguish from the incoming one
        sp = Spheroid{SoundSoft}(; equatorialRadius=0.6, polarRadius=1.0)
        md = quiet(() -> SS.modes(sp, planeWave))
        coefficients = md.A .* md.b

        for n̂ in (SVector(0.0, 0.0, 1.0), normalize(SVector(1.0, 2.0, -0.5)))
            deviation = map((true, false)) do outgoing
                Q = [SS.seriesvalue(sp, md, r * n̂, coefficients, outgoing) * r * cis(2.0 * r) for r in (10.0, 14.0, 18.0, 20.0)]
                maximum(abs.(Q .- Q[1])) / abs(Q[1])
            end

            @test deviation[1] < 1e-1
            @test deviation[2] > 3e-1
        end
    end

    @testset "Reference truncation" begin

        # the automatic settings against a deliberately generous truncation
        sp = Spheroid{SoundHard}(; equatorialRadius=0.6, polarRadius=1.0, axis=SVector(0.2, -0.3, 1.0))
        points = [r * normalize(d) for r in (1.5, 20.0) for d in (SVector(0.0, 0.0, 1.0), SVector(1.0, 2.0, -0.5))]

        reference = quiet(() -> SS.modes(sp, planeWave; M=14, N=40))
        expected = [scatteredfield(sp, planeWave, reference, p, Pressure(nothing)) for p in points]

        @test scatteredfield(sp, planeWave, Pressure(points)) ≈ expected rtol = 1e-9
    end

    @testset "Far field" begin

        # the far field is the limit of r exp(j k r) p, the deviation decreasing like 1 / r
        sp = Spheroid{SoundHard}(; equatorialRadius=0.6, polarRadius=1.0, axis=SVector(0.2, -0.3, 1.0))
        md = quiet(() -> SS.modes(sp, planeWave))

        for d in (normalize(SVector(0.3, 0.2, 1.0)), SVector(1.0, 0.0, 0.0), normalize(SVector(-0.5, 1.0, -0.4)))
            far = scatteredfield(sp, planeWave, md, FarField([d]))[1]
            deviation(r) = abs(scatteredfield(sp, planeWave, md, Pressure([r * d]))[1] * r * cis(2.0 * r) - far) / abs(far)

            @test deviation(400.0) < 1e-2
            @test deviation(800.0) / deviation(400.0) ≈ 0.5 rtol = 0.1
        end
    end

    @testset "Neumann trace at the tips" begin

        # the tips lie on the axis, where the trace is a limit; the average over two antipodal points next to it
        # approaches the axis like ρ²
        sp = Spheroid{SoundHard}(; equatorialRadius=0.6, polarRadius=1.0)
        md = quiet(() -> SS.modes(sp, planeWave))

        for tip in (1.0, -1.0)
            P = locate(sp, tip, 0.0)
            normal = outwardNormal(sp, P)

            @test normal == SVector(0.0, 0.0, tip)

            onAxis = scatteredfield(sp, planeWave, md, PressureNormalGradient([P], [normal]))[1]

            for (ρ, tolerance) in ((1e-3, 1e-5), (1e-4, 1e-7))
                η = tip * sqrt(1 - (ρ / equatorialRadius(sp))^2)
                near = [locate(sp, η, φ) for φ in (0.7, 0.7 + π)]
                average = sum(scatteredfield(sp, planeWave, md, PressureNormalGradient(near, outwardNormals(sp, near)))) / 2

                @test abs(average - onAxis) / abs(onAxis) < tolerance
            end

            @test abs(field(sp, planeWave, md, PressureNormalGradient([P], [normal]))[1]) / abs(onAxis) < 1e-10
        end
    end

    @testset "Nearby source and interior" begin

        # a monopole close to a tip requires the degree to be raised
        sp = Spheroid{SoundSoft}(; equatorialRadius=0.4, polarRadius=1.0)
        ex = SS.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.05, 0.0, 1.6))

        N₀ = ceil(Int, SS.spheroidalParameter(sp, ex) * SS.normalizedCircumradius(sp)) + 15
        md = quiet(() -> SS.modes(sp, ex))

        points = [locate(sp, η, φ) for (η, φ) in grid]

        @test md.N > N₀
        @test maximum(abs.(field(sp, ex, md, PressureTrace(points)))) / maximum(abs.(field(ex, PressureTrace(points)))) < 1e-10

        # the total pressure vanishes inside
        @test all(iszero, field(sp, ex, md, Pressure([SVector(0.0, 0.0, 0.9), SVector(0.3, 0.0, 0.0)])))
    end
end
