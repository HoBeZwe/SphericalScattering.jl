
# The modal machinery of the spheroids that is independent of the shape: the truncation, the preconditions, and
# sources close to the scatterer. It is exercised with oblate spheroids and discs.

@testitem "Automatic truncation" setup = [Setup] begin

    freq(k) = k * c / (2π)

    quantity = Pressure(nothing)

    κ = 2.0
    sp = Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)

    @testset "Measured order" begin

        # the order is measured from the azimuthal spectrum of the incident field, so that a nearly
        # axial excitation costs far fewer orders than a grazing one
        orders = map((SVector(0.0, 0.0, 1.0), normalize(SVector(0.4, 0.3, 1.0)), SVector(1.0, 0.0, 0.0))) do d
            SphericalScattering.modes(sp, SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=d)).M
        end

        @test orders[1] == 0            # an axial plane wave is rotationally symmetric
        @test orders[1] < orders[2] < orders[3]
    end

    @testset "Derived degree and overrides" begin

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

        # the degree follows from the counterpart of ka
        md = SphericalScattering.modes(sp, ex)

        @test md.N == ceil(Int, SphericalScattering.spheroidalParameter(sp, ex) * sqrt(1 + sp.ξ₀^2)) + 15
        @test size(md.A) == (2 * md.M + 1, md.N + 1)
        @test size(md.b) == size(md.A)

        # `nmax` of the parameters overrides it
        @test SphericalScattering.modes(sp, ex; parameter=Parameter(12, 1e-12)).N == 12

        # and an explicit truncation is honoured
        mdExplicit = SphericalScattering.modes(sp, ex; M=4, N=9)

        @test mdExplicit.M == 4
        @test mdExplicit.N == 9

        @test_throws ErrorException SphericalScattering.modes(sp, ex; M=12, N=4)
    end

    @testset "Accuracy of the automatic settings" begin

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

        points = [r * n̂ for r in (3.0, 20.0) for n̂ in (SVector(0.0, 0.0, 1.0), normalize(SVector(1.0, 2.0, -0.5)))]

        # a deliberately generous truncation serves as the reference
        reference = SphericalScattering.modes(sp, ex; M=8, N=40)
        expected = [scatteredfield(sp, ex, reference, point, quantity) for point in points]

        @test scatteredfield(sp, ex, Pressure(points)) ≈ expected rtol = 1e-9
    end

    @testset "Standard signatures" begin

        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))
        points = [SVector(0.0, 0.0, 3.0), SVector(2.0, 1.0, -1.5)]

        F = scatteredfield(sp, ex, Pressure(points))
        G = field(sp, ex, Pressure(points))

        @test size(F) == size(points)
        @test all(isfinite, F)
        @test G ≈ field(ex, Pressure(points)) + F rtol = 1e-12
    end

    @testset "Boundary conditions with automatic settings" begin

        # the end-to-end check of the automatic path: no truncation is given anywhere
        ex = SphericalScattering.Acoustic.planeWave(; frequency=freq(κ), direction=normalize(SVector(0.4, 0.3, 1.0)))

        locate(s, ξ, η, φ) = SphericalScattering.frame(s) * SphericalScattering.obl2cart(SVector(ξ, η, φ), s.semifocal)

        for scatterer in (Spheroid{SoundSoft}(; equatorialRadius=sqrt(1.25), polarRadius=0.5), Disc(SoundSoft; radius=1.0))
            md = SphericalScattering.modes(scatterer, ex)

            for (η, φ) in ((0.3, 0.0), (-0.5, 1.1), (0.8, 3.0))
                point = locate(scatterer, scatterer.ξ₀, η, φ)

                @test abs(field(ex, point, quantity) + scatteredfield(scatterer, ex, md, point, quantity)) /
                      abs(field(ex, point, quantity)) < 1e-10
            end
        end
    end

    @testset "Degree raised for a nearby source" begin

        # A monopole above a disc requires more degrees than the size of the disc suggests. The degree is raised
        # until the omitted modes are negligible on the surface; a given degree is used as it is, the check of
        # the truncation reporting its insufficiency.
        d = Disc(SoundSoft; radius=1.0)
        ex = SphericalScattering.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.2, 0.1, 0.6))

        N₀ = ceil(Int, SphericalScattering.spheroidalParameter(d, ex) * SphericalScattering.normalizedCircumradius(d)) + 15

        surface = [SphericalScattering.cartesianCoordinates(d, SVector(0.0, η, φ)) for η in (0.15, 0.5, 0.85) for φ in (0.3, 2.0, 4.1)]
        incident = field(ex, PressureTrace(surface))
        residual(md) = maximum(abs.(incident .+ scatteredfield(d, ex, md, PressureTrace(surface)))) / maximum(abs.(incident))

        printed(f) = mktemp() do path, io
            md = redirect_stdout(f, io)
            close(io)
            return md, read(path, String)
        end

        md, message = printed(() -> SphericalScattering.modes(d, ex))

        @test md.N > N₀
        @test residual(md) < 1e-10
        @test isempty(message)

        mdGiven, messageGiven = printed(() -> SphericalScattering.modes(d, ex; N=N₀))

        @test mdGiven.N == N₀
        @test residual(mdGiven) > 1e-7 # which the check reports
        @test occursin("truncation may be insufficient", messageGiven)
    end
end


@testitem "Spheroid interior and preconditions" setup = [Setup] begin

    freq(k) = k * c / (2π)

    sp = Spheroid{SoundHard}(; equatorialRadius=sqrt(1.25), polarRadius=0.5)

    @testset "isinside" begin

        # the predicate is answered by the scatterer in terms of its own geometry
        @test SphericalScattering.isinside(HardSphere(; radius=1.0), SVector(0.0, 0.0, 0.5))
        @test !SphericalScattering.isinside(HardSphere(; radius=1.0), SVector(0.0, 0.0, 1.5))

        # along the axis the spheroid reaches to 0.5, in the equatorial plane to sqrt(1.25)
        @test SphericalScattering.isinside(sp, SVector(0.0, 0.0, 0.2))
        @test !SphericalScattering.isinside(sp, SVector(0.0, 0.0, 0.9))
        @test SphericalScattering.isinside(sp, SVector(1.0, 0.0, 0.0))
        @test !SphericalScattering.isinside(sp, SVector(1.2, 0.0, 0.0))

        # a disc has no interior
        d = Disc(SoundHard; radius=1.0)

        @test !SphericalScattering.isinside(d, SVector(0.5, 0.0, 0.0))
        @test !SphericalScattering.isinside(d, SVector(0.0, 0.0, 0.2))
    end

    @testset "Monopole inside the scatterer" begin

        # the expansion of the incident field assumes the monopole to be outside. The check used to
        # ask the sphere for its radius, which a spheroid does not have, and was not reached at all by
        # the spheroidal series
        ex = SphericalScattering.Acoustic.monopole(; position=SVector(0.0, 0.0, 0.2), frequency=freq(2.0))

        err = ErrorException("The monopole has to be located outside the scatterer, as the expansion of its field assumes.")

        @test_throws err scatteredfield(HardSphere(; radius=1.0), ex, Pressure([SVector(0.0, 0.0, 3.0)]))
        @test_throws err scatteredfield(sp, ex, Pressure([SVector(0.0, 0.0, 3.0)]))
        @test_throws err SphericalScattering.modes(sp, ex)

        # a monopole outside is accepted
        exOutside = SphericalScattering.Acoustic.monopole(; position=SVector(0.5, -0.3, 4.0), frequency=freq(2.0))

        @test_nowarn SphericalScattering.modes(sp, exOutside)
    end
end


@testitem "Source close to a spheroid" setup = [Setup] begin

    freq(k) = k * c / (2π)
    quiet(f) = redirect_stdout(f, devnull) # the check of the truncation may report

    @testset "Source close to the scatterer" begin

        # The expansion of the field of a monopole holds between the scatterer and the monopole; the closer the
        # monopole, the more slowly it converges, so that a source closer to a disc than its radius requires
        # more degrees than the size of the disc suggests
        d = Disc(SoundSoft; radius=1.0)
        ex = SphericalScattering.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.2, 0.1, 0.3))

        md = quiet(() -> SphericalScattering.modes(d, ex; N=50))

        surface = [SphericalScattering.cartesianCoordinates(d, SVector(0.0, η, φ)) for η in (0.2, 0.5, 0.8) for φ in (0.3, 2.0)]
        incident = field(ex, PressureTrace(surface))

        @test maximum(abs.(incident .+ scatteredfield(d, ex, md, PressureTrace(surface)))) / maximum(abs.(incident)) < 1e-6

        # a plane wave has no source at a finite distance
        @test SphericalScattering.sourceCoordinate(d, SphericalScattering.Acoustic.planeWave(; frequency=freq(1.5))) == Inf
    end

    @testset "Projection between the scatterer and the source" begin

        # the projection, the fallback for excitations without an analytic expansion, has to place its surface
        # between the scatterer and the source: for a monopole closer to a disc than the default surface ξ = 0.5
        # the surface is moved, and a given surface reaching the source is rejected
        d = Disc(SoundSoft; radius=1.0)
        ex = SphericalScattering.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.2, 0.1, 0.3))

        ξs = SphericalScattering.sourceCoordinate(d, ex)

        @test ξs < SphericalScattering.projectionCoordinate(d)
        @test_throws ErrorException SphericalScattering.projectedCoefficients(d, ex, 20; ξ=1.2 * ξs)
    end
end


@testitem "Analytic and projected incident coefficients" setup = [Setup] begin

    SS = SphericalScattering

    freq(k) = k * c / (2π)

    # The analytic coefficients rest on the conventions of the spheroidal wave functions: their normalization, the
    # Condon-Shortley phase, and the normalization of the radial functions. The projection of the incident field
    # requires none of these, which makes it an independent check. Its coefficients are accurate only as far as
    # they matter on the surface of the projection, so the modes are compared by their size there.
    function deviation(sp, ex, N; quadrature...)
        ξp = min(SS.projectionCoordinate(sp), (sp.ξ₀ + SS.sourceCoordinate(sp, ex)) / 2)
        c = SS.spheroidalParameter(sp, ex)

        analytic = SS.incidentCoefficients(sp, ex, N)
        projected = SS.projectedCoefficients(sp, ex, N; quadrature...)

        weight = zeros(2N + 1, N + 1)
        for m in (-N):N
            R₁ = SS.rmn(abs(m), abs(m):N, c, [ξp]; spheroid=SS.shape(sp), kind=1).value
            for (k, n) in enumerate(abs(m):N)
                weight[m + N + 1, n + 1] = abs(R₁[1, k]) * sqrt(SS.angularNorm(abs(m), n))
            end
        end

        return maximum(abs.(analytic .- projected) .* weight) / maximum(abs.(analytic) .* weight)
    end

    scatterers = (
        Spheroid{SoundSoft}(; equatorialRadius=1.2, polarRadius=0.7, axis=SVector(0.3, -0.4, 1.0)),
        Spheroid{SoundHard}(; equatorialRadius=0.6, polarRadius=1.0, axis=SVector(0.2, -0.3, 1.0)),
        Disc(SoundSoft; radius=1.0),
    )
    excitations = (
        SS.Acoustic.planeWave(; frequency=freq(2.0), direction=normalize(SVector(0.4, 0.3, 1.0))),
        SS.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.5, -1.0, 2.5)),
    )

    for sp in scatterers, ex in excitations
        @test deviation(sp, ex, 18) < 1e-12
    end

    # --- a monopole close to a disc, the projection surface lying between the two. The field varies rapidly on
    #     that surface, so that the default quadrature of the projection limits the agreement to about 1e-10;
    #     refined, it is exact, which the analytic coefficients are regardless
    exClose = SS.Acoustic.monopole(; frequency=freq(1.5), position=SVector(0.2, 0.1, 0.3))
    N = 30

    @test deviation(Disc(SoundSoft; radius=1.0), exClose, N) < 1e-9
    @test deviation(Disc(SoundSoft; radius=1.0), exClose, N; nη=4N + 32, nφ=8N + 16) < 1e-14

    # --- the closed form of the norm, and the scattering coefficients of the batched path against those of a single mode
    sp, ex = scatterers[1], excitations[1]
    c = SS.spheroidalParameter(sp, ex)
    η, w = SS.gausslegendre(60)

    for m in 0:4
        S = SS.smn(m, m:12, c, η; spheroid=SS.shape(sp), normalize=false).value
        @test [sum(w .* view(S, :, k) .^ 2) for k in 1:(13 - m)] ≈ [SS.angularNorm(m, n) for n in m:12] rtol = 1e-13
    end

    md = SS.modes(sp, ex)
    for m in (-md.M):(md.M), n in abs(m):(md.N)
        @test md.b[m + md.M + 1, n + 1] ≈ SS.scatterCoeff(sp, ex, abs(m), n) rtol = 1e-13 # batched and single calls agree to the last digit
    end
end
