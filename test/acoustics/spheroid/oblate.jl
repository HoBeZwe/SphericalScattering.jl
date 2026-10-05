
# The oblate spheroid and its degenerate case, the disc. The machinery shared with the prolate spheroid, which
# is independent of the shape, is tested in `machinery.jl`; the coordinates in `geometry/coordinateTransforms.jl`.

@testitem "Oblate spheroid" setup = [Setup] begin

    @testset "Types" begin

        sp = Spheroid{SoundHard}(; equatorialRadius=2.0, polarRadius=1.0)

        @test sp isa OblateSpheroid{SoundHard,Float64}
        @test sp isa Scatterer{SoundHard}
        @test !(sp isa SphericalScattering.Sphere) # the spheroidal coordinates degenerate for a sphere
        @test sp.semifocal ≈ sqrt(3.0) rtol = 1e-14
        @test equatorialRadius(sp) ≈ 2.0 rtol = 1e-14
        @test polarRadius(sp) ≈ 1.0 rtol = 1e-14
        @test !isdisc(sp)
        @test sp.axis ≈ SVector(0.0, 0.0, 1.0)

        # --- the axis is normalized
        @test norm(Spheroid{SoundSoft}(; equatorialRadius=2.0, polarRadius=1.0, axis=SVector(1.0, 2.0, 3.0)).axis) ≈ 1.0

        # --- the disc is the degenerate spheroid
        d = Disc(SoundHard; radius=1.5)

        @test d isa OblateSpheroid{SoundHard,Float64}
        @test isdisc(d)
        @test d.ξ₀ == 0.0
        @test d.semifocal ≈ 1.5 rtol = 1e-14
        @test equatorialRadius(d) ≈ 1.5 rtol = 1e-14
        @test polarRadius(d) ≈ 0.0 atol = 1e-15

        # --- the spheroidal coordinates degenerate for a sphere, hence a ≠ b is required; a > b is an oblate shape,
        #     a < b a prolate one, which the oblate constructor rejects
        @test_throws ErrorException Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=1.0)
        @test Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=2.0) isa ProlateSpheroid
        @test_throws ErrorException OblateSpheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=2.0)
        @test_throws ErrorException Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=-0.5)
    end

    @testset "Spheroidal wave functions" begin

        # --- the outgoing radial function has to behave like exp(-im c ξ) / (c ξ) for large ξ, and must
        #     clearly NOT behave like exp(+im c ξ). This pins the kind of the radial function: the
        #     spheroidal literature calls the outgoing function R⁽³⁾, which belongs to the opposite time
        #     convention, so that picking that kind would yield an incoming wave.
        c = 20.0
        ξs = (8.0, 9.0, 10.0)

        outgoing = [SphericalScattering.spheroidalRadial(:oblate, 0, 1, c, ξ; outgoing=true).value for ξ in ξs]

        spread(s) = begin
            ratios = [outgoing[i] / (cis(s * c * ξs[i]) / (c * ξs[i])) for i in eachindex(ξs)]
            maximum(abs.(ratios .- ratios[1])) / abs(ratios[1])
        end

        @test spread(-1) < 1e-1   # matches the outgoing wave
        @test spread(+1) > 1.0    # and is far from the incoming one

        # --- the Wronskian of the two real solutions is known: R₁R₂' - R₂R₁' = 1 / (c (ξ² + 1)).
        #     Reconstructing R₂ from the regular and the outgoing function tests the adapter itself.
        for c in (0.5, 5.0, 100.0), m in (0, 2), n in (m, m + 3), ξ in (0.0, 0.4, 2.0)
            R₁ = SphericalScattering.spheroidalRadial(:oblate, m, n, c, ξ)
            Rₒ = SphericalScattering.spheroidalRadial(:oblate, m, n, c, ξ; outgoing=true)

            R₂ = im * (Rₒ.value - R₁.value)              # since Rₒ = R₁ - im R₂
            dR₂ = im * (Rₒ.derivative - R₁.derivative)

            @test real(R₁.value * dR₂ - R₂ * R₁.derivative) ≈ 1 / (c * (ξ^2 + 1)) rtol = 1e-10
        end

        # --- the angular function reduces to the Legendre function as c → 0
        for m in (0, 1), n in (m, m + 2), η in (-0.6, 0.3)
            @test SphericalScattering.spheroidalAngular(:oblate, m, n, 1e-7, η).value ≈ SphericalScattering.Plm(η, n, m) rtol = 1e-6
        end
    end

    @testset "Scattering coefficients" begin

        f = 1e8
        ex = SphericalScattering.Acoustic.planeWave(; frequency=f)

        sp = Spheroid{SoundHard}(; equatorialRadius=2.0, polarRadius=1.0)
        sps = Spheroid{SoundSoft}(; equatorialRadius=2.0, polarRadius=1.0)

        # --- for a closed spheroid every mode scatters
        for m in (0, 1), n in (m, m + 1, m + 2)
            @test isfinite(SphericalScattering.scatterCoeff(sp, ex, m, n))
            @test isfinite(SphericalScattering.scatterCoeff(sps, ex, m, n))
            @test !iszero(SphericalScattering.scatterCoeff(sp, ex, m, n))
            @test !iszero(SphericalScattering.scatterCoeff(sps, ex, m, n))
        end

        # --- on a disc the radial functions split by parity: the Neumann problem is carried by the modes
        #     with odd n - m, the Dirichlet problem by those with even n - m
        dHard = Disc(SoundHard; radius=1.0)
        dSoft = Disc(SoundSoft; radius=1.0)

        for m in (0, 1), n in m:(m + 3)
            bHard = SphericalScattering.scatterCoeff(dHard, ex, m, n)
            bSoft = SphericalScattering.scatterCoeff(dSoft, ex, m, n)

            if iseven(n - m)
                @test iszero(bHard)
                @test !iszero(bSoft)
            else
                @test !iszero(bHard)
                @test iszero(bSoft)
            end
        end

        # --- the coefficients do not depend on the excitation
        exMono = SphericalScattering.Acoustic.monopole(; position=SVector(0.0, 0.0, 5.0), frequency=f)

        @test SphericalScattering.scatterCoeff(sp, exMono, 0, 1) === SphericalScattering.scatterCoeff(sp, ex, 0, 1)
    end
end


@testitem "Oblate expansion" setup = [Setup] begin

    freq(k) = k * c / (2π)   # the frequency belonging to a desired wavenumber

    quantity = Pressure(nothing)

    # a point of the surface ξ = const, in Cartesian coordinates of the global frame
    locate(sp, ξ, η, φ) = SphericalScattering.frame(sp) * SphericalScattering.obl2cart(SVector(ξ, η, φ), sp.semifocal)

    surface = [(0.3, 0.0), (-0.5, 1.1), (0.8, 3.0), (0.0, 2.0), (-0.95, 4.5)]

    @testset "Frame of the spheroid" begin

        for axis in (SVector(0.0, 0.0, 1.0), normalize(SVector(1.0, 1.0, 1.0)), SVector(0.0, 0.0, -1.0))
            sp = Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5, axis=axis)
            R = SphericalScattering.frame(sp)

            @test R * R' ≈ I atol = 1e-14          # orthogonal
            @test det(R) ≈ 1.0 atol = 1e-14        # and a rotation
            @test R * SVector(0.0, 0.0, 1.0) ≈ axis atol = 1e-14
        end
    end

    @testset "Incident re-expansion" begin

        sp = Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)   # semifocal 1, ξ₀ = 0.5
        points = [SVector(0.3, -0.4, 1.6), SVector(0.0, 0.0, 2.5), SVector(1.1, 0.0, 0.0)]

        # --- projecting the incident field onto the angular functions and reconstructing it elsewhere
        #     converges with the truncation: this validates the expansion without any closed-form
        #     coefficients, and hence independently of the normalization of the angular functions
        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))

        errors = map(
            ((M, N),) -> begin
                md = SphericalScattering.modes(sp, ex; M=M, N=N, ξ=3.0)
                maximum(
                    abs(SphericalScattering.seriesvalue(sp, md, point, md.A, false) - field(ex, point, quantity)) /
                    abs(field(ex, point, quantity)) for point in points
                )
            end,
            ((4, 10), (8, 20)),
        )

        @test errors[2] < errors[1] / 100   # the error drops steeply with the truncation
        @test errors[2] < 1e-6

        # --- a monopole is expanded just as well, the projection surface lying inside the source
        exMono = SphericalScattering.Acoustic.monopole(; position=SVector(0.5, -0.3, 4.0), frequency=freq(2.0))
        md = SphericalScattering.modes(sp, exMono; M=8, N=20)

        for point in (SVector(1.1, 0.0, 0.0), SVector(0.3, -0.4, 0.9))
            @test SphericalScattering.seriesvalue(sp, md, point, md.A, false) ≈ field(exMono, point, quantity) rtol = 1e-6
        end

        # --- an arbitrarily oriented spheroid is handled by transforming into its frame
        spTilted = Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5, axis=normalize(SVector(1.0, 1.0, 1.0)))
        mdTilted = SphericalScattering.modes(spTilted, ex; M=8, N=20, ξ=3.0)

        for point in points
            @test SphericalScattering.seriesvalue(spTilted, mdTilted, point, mdTilted.A, false) ≈ field(ex, point, quantity) rtol =
                1e-5
        end
    end

    @testset "Boundary conditions" begin

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))

        # --- the total pressure vanishes on a sound-soft spheroid. The incident field is known in closed
        #     form, so this also verifies the coefficients of the incident expansion
        spSoft = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
        mdSoft = SphericalScattering.modes(spSoft, ex; M=8, N=20)

        for (η, φ) in surface
            point = locate(spSoft, spSoft.ξ₀, η, φ)

            @test abs(field(ex, point, quantity) + scatteredfield(spSoft, ex, mdSoft, point, Pressure(nothing))) /
                  abs(field(ex, point, quantity)) < 1e-7
        end

        # --- the radial derivative of the total pressure vanishes on a sound-hard spheroid
        spHard = Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
        mdHard = SphericalScattering.modes(spHard, ex; M=8, N=20)

        h = 1e-5

        for (η, φ) in surface
            total(ξ) = begin
                point = locate(spHard, ξ, η, φ)
                field(ex, point, quantity) + scatteredfield(spHard, ex, mdHard, point, Pressure(nothing))
            end

            derivative = (-3 * total(spHard.ξ₀) + 4 * total(spHard.ξ₀ + h) - total(spHard.ξ₀ + 2 * h)) / (2 * h)
            scale = abs((total(spHard.ξ₀ + h) - total(spHard.ξ₀)) / h) + abs(total(spHard.ξ₀))

            @test abs(derivative) / scale < 1e-6
        end

        # --- the disc is the degenerate spheroid ξ₀ = 0 and needs no special treatment
        dSoft = Disc(SoundSoft; radius=1.0)
        mdDisc = SphericalScattering.modes(dSoft, ex; M=8, N=20)

        for (η, φ) in surface
            point = locate(dSoft, 0.0, η, φ)

            @test abs(field(ex, point, quantity) + scatteredfield(dSoft, ex, mdDisc, point, Pressure(nothing))) /
                  abs(field(ex, point, quantity)) < 1e-7
        end

        # --- and a monopole at an arbitrary position works as well
        exMono = SphericalScattering.Acoustic.monopole(; position=SVector(0.5, -0.3, 4.0), frequency=freq(2.0))
        mdMono = SphericalScattering.modes(spSoft, exMono; M=8, N=20)

        for (η, φ) in surface
            point = locate(spSoft, spSoft.ξ₀, η, φ)

            @test abs(field(exMono, point, quantity) + scatteredfield(spSoft, exMono, mdMono, point, Pressure(nothing))) /
                  abs(field(exMono, point, quantity)) < 1e-9
        end
    end

    @testset "Array interface" begin

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))
        sp = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
        md = SphericalScattering.modes(sp, ex; M=6, N=14)

        points = [locate(sp, 1.5, η, φ) for (η, φ) in surface]

        F = scatteredfield(sp, ex, md, Pressure(points))
        G = field(sp, ex, md, Pressure(points))

        @test size(F) == size(points)
        @test all(isfinite, F)
        @test F[1] ≈ scatteredfield(sp, ex, md, points[1], Pressure(nothing)) rtol = 1e-14
        @test G ≈ field(ex, Pressure(points)) + F rtol = 1e-14
    end
end


@testitem "Oblate independent validation" setup = [Setup] begin

    freq(k) = k * c / (2π)

    quantity = Pressure(nothing)

    # The boundary conditions cannot validate the scattered field on their own: the scattering
    # coefficients are defined by them, so that the residual reduces to the error of the incident
    # re-expansion on the surface. In particular, an incoming wave would fulfill them just as well.
    # The two tests below close that gap.

    @testset "Sphere bridge" begin

        # As the eccentricity vanishes the oblate spheroid degenerates into a sphere: the semifocal
        # distance tends to zero and ξ₀ to infinity, with c ξ₀ → k a. Since the geometries then differ
        # by the eccentricity, so do the fields, which ties the spheroidal series to the independently
        # validated series of the sphere.
        a = 1.0
        κ = 2.0

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=SVector(0.0, 0.0, 1.0))

        points = [SVector(0.0, 0.0, 3.0), SVector(2.0, 1.0, -1.5), SVector(-3.0, 0.0, 0.0)]

        for (BC, SphereType) in ((SoundSoft, SoftSphere), (SoundHard, HardSphere))

            reference = [scatteredfield(SphereType(; radius=a), ex, Pressure([point]))[1] for point in points]

            deviation = map((1e-2, 1e-4)) do δ
                sp = Spheroid{BC}(; equatorialRadius=a * (1 + δ), polarRadius=a)
                md = SphericalScattering.modes(sp, ex; M=2, N=16)

                # the parametrization has to preserve c ξ₀ = κ a
                @test SphericalScattering.spheroidalParameter(sp, ex) * sp.ξ₀ ≈ κ * a rtol = 1e-12

                got = [scatteredfield(sp, ex, md, point, quantity) for point in points]

                maximum(abs.(got .- reference) ./ abs.(reference))
            end

            # the deviation is of the order of the eccentricity, hence drops by two decades
            @test deviation[2] < deviation[1] / 50
            @test deviation[2] < 1e-3
        end
    end

    @testset "Radiation condition" begin

        # The scattered field has to decay as exp(-im κ r) / r, so that p_s r exp(im κ r) approaches a
        # constant. Evaluating the series with the regular instead of the outgoing radial function is
        # rejected by a wide margin, which no boundary condition would do.
        κ = 2.0

        sp = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

        # the degree bounds the scattered field at every radius, so the automatic truncation suffices
        md = SphericalScattering.modes(sp, ex)

        coefficients = md.A .* md.b
        radii = (10.0, 14.0, 18.0, 20.0)

        for n̂ in (SVector(0.0, 0.0, 1.0), normalize(SVector(1.0, 2.0, -0.5)))

            deviation = map((true, false)) do outgoing
                Q = [SphericalScattering.seriesvalue(sp, md, r * n̂, coefficients, outgoing) * r * cis(κ * r) for r in radii]

                maximum(abs.(Q .- Q[1])) / abs(Q[1])
            end

            @test deviation[1] < 1e-1   # the outgoing solution leaves p_s r exp(im κ r) nearly constant
            @test deviation[2] > 3e-1   # the regular one does not
        end
    end
end


@testitem "Oblate traces" setup = [Setup] begin

    freq(k) = k * c / (2π)

    quantity = Pressure(nothing)
    trace = PressureTrace(nothing)

    κ = 2.0
    ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

    locate(s, ξ, η, φ) = SphericalScattering.frame(s) * SphericalScattering.obl2cart(SVector(ξ, η, φ), s.semifocal)

    # η = 0 is the equator of a spheroid, but the rim of a disc, where the basis degenerates
    surface = [(0.3, 0.0), (-0.5, 1.1), (0.8, 3.0), (0.0, 2.0)]
    surfaceOffRim = [(0.3, 0.0), (-0.5, 1.1), (0.8, 3.0)]

    sp = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
    md = SphericalScattering.modes(sp, ex)

    @testset "Outward normal" begin

        # the outward normal of a spheroid is ê_ξ, which is perpendicular to both surface directions
        for (η, φ) in surface
            point = locate(sp, sp.ξ₀, η, φ)
            n̂ = outwardNormal(sp, point)

            @test norm(n̂) ≈ 1.0 rtol = 1e-14

            ~, êη, êφ = SphericalScattering.oblateBasis(SVector(sp.ξ₀, η, φ), sp.semifocal)
            R = SphericalScattering.frame(sp)

            @test abs(dot(n̂, R * êη)) < 1e-13
            @test abs(dot(n̂, R * êφ)) < 1e-13
        end

        # it is not parallel to the position vector, except at the equator and at the poles, where the
        # two coincide by symmetry
        for (η, φ) in ((0.3, 0.0), (-0.5, 1.1), (0.8, 3.0))
            point = locate(sp, sp.ξ₀, η, φ)

            @test abs(dot(outwardNormal(sp, point), normalize(point))) < 0.999
        end

        for (η, φ) in ((0.0, 2.0),)
            point = locate(sp, sp.ξ₀, η, φ)

            @test abs(dot(outwardNormal(sp, point), normalize(point))) ≈ 1.0 rtol = 1e-12
        end

        # on a disc it is the face normal
        d = Disc(SoundHard; radius=1.0)

        for (η, φ) in surfaceOffRim
            @test outwardNormal(d, locate(d, 0.0, η, φ)) ≈ SVector(0.0, 0.0, 1.0) atol = 1e-14
        end

        # the convenience constructor fills them in
        points = [locate(sp, sp.ξ₀, η, φ) for (η, φ) in surface]

        @test PressureNormalGradient(sp, points).normals ≈ outwardNormals(sp, points)
        @test !(PressureNormalGradient(sp, points).normals ≈ PressureNormalGradient(points).normals)
    end

    @testset "Dirichlet trace" begin

        # the trace is the scattered pressure evaluated on the surface
        for (η, φ) in surface
            point = locate(sp, sp.ξ₀, η, φ)

            @test scatteredfield(sp, ex, md, point, trace) ≈ scatteredfield(sp, ex, md, point, quantity) rtol = 1e-12
        end

        # only the direction of the location matters
        for (η, φ) in surface
            onSurface = locate(sp, sp.ξ₀, η, φ)

            @test scatteredfield(sp, ex, md, locate(sp, 3.0, η, φ), trace) ≈ scatteredfield(sp, ex, md, onSurface, trace) rtol = 1e-12
        end
    end

    @testset "Gradient" begin

        # The partial derivatives are compared with differences taken along the surface, so that
        # nothing is sampled in the interior: the one with respect to ξ is one-sided and outward, the
        # tangential ones are central and stay on the surface.
        h = 1e-6
        coefficients = md.A .* md.b
        R = SphericalScattering.frame(sp)

        for (η, φ) in surfaceOffRim

            point = locate(sp, sp.ξ₀, η, φ)

            gradient = SphericalScattering.surfaceseries(sp, md, point, coefficients).gradient
            metric = SphericalScattering.oblateMetric(SVector(sp.ξ₀, η, φ), sp.semifocal)
            êξ, êη, êφ = SphericalScattering.oblateBasis(SVector(sp.ξ₀, η, φ), sp.semifocal)

            pressure(ξ) = SphericalScattering.seriesvalue(sp, md, locate(sp, ξ, η, φ), coefficients, true)
            surfaceValue(ηη, φφ) = SphericalScattering.surfaceseries(sp, md, locate(sp, sp.ξ₀, ηη, φφ), coefficients).value

            ∂ξ = (-3 * pressure(sp.ξ₀) + 4 * pressure(sp.ξ₀ + h) - pressure(sp.ξ₀ + 2 * h)) / (2 * h)
            ∂η = (surfaceValue(η + h, φ) - surfaceValue(η - h, φ)) / (2 * h)
            ∂φ = (surfaceValue(η, φ + h) - surfaceValue(η, φ - h)) / (2 * h)

            @test dot(R * êξ, gradient) * metric[1] ≈ ∂ξ rtol = 1e-7
            @test dot(R * êη, gradient) * metric[2] ≈ ∂η rtol = 1e-7
            @test dot(R * êφ, gradient) * metric[3] ≈ ∂φ rtol = 1e-7
        end
    end

    @testset "Boundary conditions" begin

        for (BC, scatterers) in (
            (SoundSoft, (Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5), Disc(SoundSoft; radius=1.0))),
            (SoundHard, (Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5), Disc(SoundHard; radius=1.0))),
        )
            for scatterer in scatterers

                points = [locate(scatterer, scatterer.ξ₀, η, φ) for (η, φ) in surfaceOffRim]
                m = SphericalScattering.modes(scatterer, ex)

                if BC === SoundSoft
                    # the total pressure vanishes on the surface
                    total = field(ex, PressureTrace(points)) .+ scatteredfield(scatterer, ex, m, PressureTrace(points))
                    reference = abs.(field(ex, PressureTrace(points)))
                else
                    # the normal derivative of the total pressure vanishes, the outward normal being ê_ξ
                    q = PressureNormalGradient(scatterer, points)
                    total = field(ex, q) .+ scatteredfield(scatterer, ex, m, q)
                    reference = abs.(field(ex, q))
                end

                @test maximum(abs.(total) ./ reference) < 1e-8
            end
        end
    end

    @testset "Disc faces" begin

        # the face cannot be recovered from Cartesian coordinates, but the parity of the angular
        # functions relates the two: γ₀ is even in η for a sound-soft and odd for a sound-hard disc
        for (BC, parity) in ((SoundSoft, +1), (SoundHard, -1))

            d = Disc(BC; radius=1.0)
            m = SphericalScattering.modes(d, ex)
            C = m.A .* m.b

            evaluate(η, φ) = sum(
                C[mm + m.M + 1, n + 1] *
                SphericalScattering.rmn(abs(mm), abs(mm):(m.N), m.c, [0.0]; spheroid=:oblate, kind=4).value[1, n - abs(mm) + 1] *
                SphericalScattering.smn(abs(mm), abs(mm):(m.N), m.c, [η]; spheroid=:oblate, normalize=false).value[
                    1, n - abs(mm) + 1
                ] *
                cis(mm * φ) for mm in (-m.M):(m.M) for n in abs(mm):(m.N)
            )

            for (η, φ) in ((0.4, 0.7), (0.8, 2.3))
                @test evaluate(-η, φ) ≈ parity * evaluate(η, φ) rtol = 1e-10
            end
        end
    end

    @testset "Array interface and total traces" begin

        points = [locate(sp, sp.ξ₀, η, φ) for (η, φ) in surface]

        for q in (PressureTrace(points), PressureNormalGradient(sp, points))
            F = scatteredfield(sp, ex, md, q)
            G = field(sp, ex, md, q)

            @test size(F) == size(points)
            @test all(isfinite, F)
            @test G ≈ field(ex, q) + F rtol = 1e-12
        end

        # and the signatures that determine the coefficients themselves
        @test scatteredfield(sp, ex, PressureTrace(points)) ≈ scatteredfield(sp, ex, md, PressureTrace(points)) rtol = 1e-9
    end
end


@testitem "Oblate far field" setup = [Setup] begin

    freq(k) = k * c / (2π)

    quantity = Pressure(nothing)
    farfield = FarField(nothing)

    κ = 2.0
    ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

    sp = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)
    md = SphericalScattering.modes(sp, ex)

    directions = [SVector(0.0, 0.0, 1.0), normalize(SVector(1.0, 2.0, -0.5)), SVector(1.0, 0.0, 0.0)]

    @testset "Limit of the pressure" begin

        # the far field omits exp(-im κ r) / r, so that the pressure has to approach it with an error
        # of the order 1 / (κ r)
        for n̂ in directions
            FF = scatteredfield(sp, ex, md, n̂, farfield)

            errors = map((100.0, 1000.0)) do r
                abs(scatteredfield(sp, ex, md, r * n̂, quantity) - FF * cis(-κ * r) / r) / abs(FF / r)
            end

            @test errors[2] < errors[1] / 5   # the error drops by a decade per decade in r
            @test errors[2] < 1e-2
        end
    end

    @testset "Determined by the direction alone" begin

        n̂ = normalize(SVector(1.0, 2.0, -0.5))
        FF = scatteredfield(sp, ex, md, n̂, farfield)

        for scale in (0.1, 1.0, 1e5)
            @test scatteredfield(sp, ex, md, scale * n̂, farfield) ≈ FF rtol = 1e-14
        end
    end

    @testset "Sphere bridge" begin

        # the far field degenerates onto that of the sphere, which ties down the factor im^(n+1) / k
        a = 1.0
        exAxial = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=SVector(0.0, 0.0, 1.0))

        for (BC, SphereType) in ((SoundSoft, SoftSphere), (SoundHard, HardSphere))

            reference = [scatteredfield(SphereType(; radius=a), exAxial, FarField([n̂]))[1] for n̂ in directions]

            deviation = map((1e-2, 1e-4)) do δ
                s = Spheroid{BC}(; equatorialRadius=a * (1 + δ), polarRadius=a)
                m = SphericalScattering.modes(s, exAxial; M=2, N=16)

                maximum(abs.([scatteredfield(s, exAxial, m, n̂, farfield) for n̂ in directions] .- reference) ./ abs.(reference))
            end

            @test deviation[2] < deviation[1] / 50
            @test deviation[2] < 1e-3
        end
    end

    @testset "Disc" begin

        for BC in (SoundSoft, SoundHard)
            d = Disc(BC; radius=1.0)
            m = SphericalScattering.modes(d, ex)

            n̂ = SVector(0.0, 0.0, 1.0)
            r = 500.0

            FF = scatteredfield(d, ex, m, n̂, farfield)

            @test isfinite(FF)
            @test scatteredfield(d, ex, m, r * n̂, quantity) ≈ FF * cis(-κ * r) / r rtol = 1e-2
        end
    end

    @testset "Array interface" begin

        F = scatteredfield(sp, ex, md, FarField(directions))

        @test size(F) == size(directions)
        @test all(isfinite, F)
        @test F ≈ [scatteredfield(sp, ex, md, n̂, farfield) for n̂ in directions] rtol = 1e-14

        # and the signature determining the coefficients itself
        @test scatteredfield(sp, ex, FarField(directions)) ≈ F rtol = 1e-9
    end

    @testset "Interior of the total pressure" begin

        # the total pressure vanishes inside, as for the spherical scatterers
        inside = [
            SphericalScattering.frame(sp) * SphericalScattering.obl2cart(SVector(0.5 * sp.ξ₀, η, φ), sp.semifocal) for
            (η, φ) in ((0.3, 0.0), (-0.5, 1.1))
        ]

        @test norm(scatteredfield(sp, ex, md, Pressure(inside))) == 0.0
        @test norm(field(sp, ex, md, Pressure(inside))) == 0.0
    end
end


@testitem "Pressure jump across a disc" setup = [Setup] begin

    freq(k) = k * c / (2π)

    trace = PressureTrace(nothing)
    jump = PressureJump(nothing)

    κ = 2.0
    ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

    locate(d, η, φ) = SphericalScattering.frame(d) * SphericalScattering.obl2cart(SVector(0.0, η, φ), d.semifocal)

    # the series evaluated explicitly at a given η, bypassing the transform, which cannot tell the faces apart
    function faceValue(md, η, φ)
        coefficients = md.A .* md.b

        return sum(
            coefficients[m + md.M + 1, n + 1] *
            SphericalScattering.rmn(abs(m), abs(m):(md.N), md.c, [0.0]; spheroid=:oblate, kind=4).value[1, n - abs(m) + 1] *
            SphericalScattering.smn(abs(m), abs(m):(md.N), md.c, [η]; spheroid=:oblate, normalize=false).value[1, n - abs(m) + 1] *
            cis(m * φ) for m in (-md.M):(md.M) for n in abs(m):(md.N)
        )
    end

    surface = [(0.4, 0.7), (0.8, 2.3), (0.6, 1.0)]   # η > 0, the face the convention refers to

    @testset "Sound-hard disc" begin

        d = Disc(SoundHard; radius=1.0)
        md = SphericalScattering.modes(d, ex)

        for (η, φ) in surface
            point = locate(d, η, φ)

            # the jump is the difference of the two faces
            @test scatteredfield(d, ex, md, point, jump) ≈ faceValue(md, η, φ) - faceValue(md, -η, φ) rtol = 1e-12

            # the surface retains the modes of odd n - m alone, whose angular functions are odd in η,
            # so that the jump is twice the trace
            @test scatteredfield(d, ex, md, point, jump) ≈ 2 * scatteredfield(d, ex, md, point, trace) rtol = 1e-12
        end
    end

    @testset "Sound-soft disc" begin

        d = Disc(SoundSoft; radius=1.0)
        md = SphericalScattering.modes(d, ex)

        # only the modes of even n - m contribute, which are the same on both faces: the pressure is
        # continuous and its jump vanishes. The unknown of a sound-soft disc is the jump of the normal
        # derivative instead
        for (η, φ) in surface
            @test scatteredfield(d, ex, md, locate(d, η, φ), jump) == 0.0
        end
    end

    @testset "Total jump and the incident field" begin

        d = Disc(SoundHard; radius=1.0)
        md = SphericalScattering.modes(d, ex)

        points = [locate(d, η, φ) for (η, φ) in surface]

        # the incident field is regular across the disc, hence the total jump is the scattered one. The
        # comparison allows for the last digit: two evaluations on several threads need not agree bitwise
        @test norm(field(ex, PressureJump(points))) == 0.0
        @test field(d, ex, PressureJump(points)) ≈ scatteredfield(d, ex, md, PressureJump(points)) rtol = 1e-14
    end

    @testset "Convention for the faces" begin

        # the faces share their Cartesian coordinates, so the jump is always reported for the face
        # whose normal is the axis: a location with η < 0 yields the same value
        d = Disc(SoundHard; radius=1.0)
        md = SphericalScattering.modes(d, ex)

        for (η, φ) in surface
            @test scatteredfield(d, ex, md, locate(d, -η, φ), jump) ≈ scatteredfield(d, ex, md, locate(d, η, φ), jump) rtol = 1e-12
        end
    end

    @testset "Array interface and closed surfaces" begin

        d = Disc(SoundHard; radius=1.0)
        md = SphericalScattering.modes(d, ex)

        points = [locate(d, η, φ) for (η, φ) in surface]

        F = scatteredfield(d, ex, md, PressureJump(points))

        @test size(F) == size(points)
        @test all(isfinite, F)
        @test scatteredfield(d, ex, PressureJump(points)) ≈ F rtol = 1e-9

        # a jump is not defined across a closed surface
        err = ErrorException("The jump of the pressure is defined across an open surface, that is, across a disc.")

        @test_throws err scatteredfield(Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5), ex, PressureJump(points))
        @test_throws err scatteredfield(HardSphere(; radius=1.0), ex, PressureJump(points))
    end
end


@testitem "Oblate axis and disc faces" setup = [Setup] begin

    freq(k) = k * c / (2π)
    quiet(f) = redirect_stdout(f, devnull) # the check of the truncation may report

    @testset "Neumann trace on the axis" begin

        # On the axis the basis vectors ê_η and ê_φ are not determined, but the surface is smooth and the
        # trace finite. Its limit is checked against the average over two antipodal points at a distance ρ
        # from the axis, which approaches the axis like ρ², the error of first order cancelling.
        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))
        locate(sp, η, φ) = SphericalScattering.frame(sp) * SphericalScattering.cartesianCoordinates(sp, SVector(sp.ξ₀, η, φ))

        arbitrary = normalize(SVector(0.3, -0.8, 0.5))

        for (sp, poles) in (
            (Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5), (1.0, -1.0)),
            (Spheroid{SoundSoft}(; equatorialRadius=1.2, polarRadius=0.7, axis=SVector(0.3, -0.4, 1.0)), (1.0, -1.0)),
            (Disc(SoundHard; radius=1.0), (1.0,)), # the faces of a disc share the center, the upper one is reported
        )
            md = quiet(() -> SphericalScattering.modes(sp, ex))
            a = equatorialRadius(sp)

            for pole in poles
                P = locate(sp, pole, 0.0)
                normal = outwardNormal(sp, P)

                for n̂ in (normal, arbitrary)
                    onAxis = scatteredfield(sp, ex, md, PressureNormalGradient([P], [n̂]))[1]

                    @test isfinite(onAxis)

                    for (ρ, tolerance) in ((1e-3, 1e-5), (1e-4, 1e-7))
                        η = pole * sqrt(1 - (ρ / a)^2)
                        near = [locate(sp, η, φ) for φ in (0.7, 0.7 + π)]
                        normals = n̂ === normal ? outwardNormals(sp, near) : [n̂, n̂]
                        average = sum(scatteredfield(sp, ex, md, PressureNormalGradient(near, normals))) / 2

                        @test abs(average - onAxis) / abs(onAxis) < tolerance
                    end
                end

                # the boundary condition holds on the axis as well
                if sp.boundary isa SoundHard
                    @test abs(field(sp, ex, md, PressureNormalGradient([P], [normal]))[1]) /
                          abs(field(ex, PressureNormalGradient([P], [normal]))[1]) < 1e-8
                end
            end
        end
    end

    @testset "Faces of a disc" begin

        # in the plane of a disc the face cannot be recovered from a location, and the upper one is reported;
        # this must not depend on the sign of a vanishing z, an artifact of rounding
        for x in (0.3, -0.3), y in (0.2, -0.2)
            @test SphericalScattering.cart2obl(SVector(x, y, -0.0), 1.0)[2] > 0
            @test SphericalScattering.cart2obl(SVector(x, y, -0.0), 1.0) == SphericalScattering.cart2obl(SVector(x, y, 0.0), 1.0)
        end

        d = Disc(SoundHard; radius=1.0)
        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0)))
        md = quiet(() -> SphericalScattering.modes(d, ex))

        points = [SVector(-0.3, -0.2, z) for z in (0.0, -0.0)]

        # the trace of a sound-hard disc is odd in η: the other face would flip its sign
        @test scatteredfield(d, ex, md, PressureTrace(points[1:1]))[1] ≈ scatteredfield(d, ex, md, PressureTrace(points[2:2]))[1] rtol =
            1e-12
    end
end
