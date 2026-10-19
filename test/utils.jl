
@testitem "Plotting points" setup = [Setup] begin

    # just check whether function calls returns no errors
    @test_nowarn sphericalGridPoints()
    @test_nowarn phiCutPoints(2)
    @test_nowarn thetaCutPoints(20)

    # --- points on the surface of a sphere, the last column repeating the first to close the surface
    pts = surfacePoints(HardSphere(; radius=2.0); resolution=10)

    @test size(pts) == (19, 37)
    @test all(p -> norm(p) ≈ 2.0, pts)
    @test pts[:, 1] ≈ pts[:, end]

    # --- on spheroids, the first row being the pole in the direction of the axis
    tilt = normalize(SVector(0.3, -0.5, 0.8))

    for sp in (
        Spheroid{SoundHard}(; equatorialRadius=1.0, polarRadius=0.4, axis=tilt),
        Spheroid{SoundSoft}(; equatorialRadius=0.5, polarRadius=1.5),
    )
        pts = surfacePoints(sp; resolution=10)
        R = SphericalScattering.frame(sp)

        @test size(pts) == (19, 37)
        @test all(p -> SphericalScattering.spheroidalCoordinates(sp, R' * p)[1] ≈ sp.ξ₀, pts)
        @test dot(pts[1, 1], sp.axis) ≈ polarRadius(sp)
        @test pts[:, 1] ≈ pts[:, end]
    end

    # --- on a disc, the rows of the two faces mirroring each other, the rim being omitted
    d = Disc(SoundHard; radius=1.0, axis=tilt)
    pts = surfacePoints(d; resolution=10)
    R = SphericalScattering.frame(d)
    coordinates = [SphericalScattering.spheroidalCoordinates(d, R' * p) for p in pts]
    η = getindex.(coordinates, 2)

    @test size(pts) == (20, 37)
    @test all(c -> c[1] < 1e-6, coordinates)
    @test all(>(0), η[1:10, :]) && all(<(0), η[11:20, :])
    @test η[1:10, :] ≈ -η[20:-1:11, :]
    @test minimum(abs, η) > 0.05 # the rim, where the Neumann trace is singular
end

@testitem "Surface densities on a disc" setup = [Setup] begin

    # the points of `surfacePoints` determine the face of a disc, also for a tilted one: the scattered pressure of a
    # sound-hard disc is odd with respect to its plane, and its total Neumann trace vanishes on both faces
    ex = Acoustic.planeWave(; frequency=0.25, embedding=Medium(1.0, 1.0), direction=normalize(SVector(1.0, 0.2, 1.0)))

    for axis in (SVector(0.0, 0.0, 1.0), normalize(SVector(0.3, -0.5, 0.8)))
        d = Disc(SoundHard; radius=1.0, axis=axis)
        pts = surfacePoints(d; resolution=15)
        n = size(pts, 1) ÷ 2

        γ₀ = scatteredfield(d, ex, PressureTrace(pts))
        @test maximum(abs.(γ₀[1:n, :] + γ₀[end:-1:(n + 1), :])) < 1e-12 * maximum(abs.(γ₀))

        q = PressureNormalGradient(d, pts)
        @test maximum(abs.(field(d, ex, q))) < 1e-6 * maximum(abs.(field(ex, q)))
    end
end

@testitem "Plotting extensions" setup = [Setup] begin

    using PlotlyJS

    # --- excite PEC sphere by Hertzian dipole
    orient = normalize(SVector(0.0, 1.0, 1.0))
    ex = HertzianDipole(; frequency=1e8, orientation=orient, position=2 * orient)
    sp = PECSphere(; radius=1.0)

    # --- evaluate fields at spherical grid points
    points_cart, points_sph = sphericalGridPoints()
    FF = scatteredfield(sp, ex, FarField(points_cart))

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=true, type="abs")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=true, type="abs")

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=true, type="theta")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=true, type="phi")

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=false, type="theta")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=false, type="phi")

    # --- evaluate fields at φ = 5° cut
    points_cart, points_sph = phiCutPoints(5) # analogously, thetaCutPoints can be used
    FF = scatteredfield(sp, ex, FarField(points_cart))

    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="polar")
    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="linear", format="polar")

    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="rectangular")
    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="linear", format="rectangular")

    # --- surface densities of acoustic scatterers
    ex = Acoustic.planeWave(; frequency=0.25, embedding=Medium(1.0, 1.0), direction=normalize(SVector(1.0, 0.0, 1.0)))
    sphere = HardSphere(; radius=1.0)

    pts = surfacePoints(sphere; resolution=30)
    γ₀ = field(sphere, ex, PressureTrace(pts))

    for kw in ((;), (type="real",), (type="imag",), (type="phase",), (scale="log", normalize=true), (clip=0.9,))
        @test length(plotsurface(sphere, γ₀, pts; kw...).plot.data) == 1
    end

    @test_throws ErrorException plotsurface(sphere, γ₀[1:2, :], pts)
    @test_throws ErrorException plotsurface(sphere, γ₀, pts; type="angle")
    @test_throws ErrorException plotsurface(sphere, γ₀, pts; type="real", scale="log")
    @test_throws ErrorException plotsurface(sphere, γ₀, pts; clip=0.0)

    oblate = Spheroid{SoundSoft}(; equatorialRadius=1.0, polarRadius=0.5)
    @test length(plotsurface(oblate, ex, PressureNormalGradient; resolution=30).plot.data) == 1

    p = plotsurface(Disc(SoundHard; radius=1.0), ex, PressureJump; resolution=30)
    @test length(p.plot.data) == 2 # the two faces of the disc
    @test p.plot.data[1][:colorbar][:title] == "|[p]|"

    p = plotsurface(sphere, ex, PressureTrace; total=false, type="real", resolution=30)
    @test p.plot.data[1][:colorbar][:title] == "Re γ₀pˢ"
end

@testitem "Radar cross section" setup = [Setup] begin

    # --- monostatic RCS
    @test_nowarn rcs(PECSphere(; radius=2.0), planeWave(; frequency=1e8))

    @test_throws ErrorException("The (monostatic) RCS is only defined for a plane-wave excitation (so far).") rcs(
        PECSphere(; radius=1.0), SphericalModeTE(; frequency=f)
    )


    # --- bistatic RCS
    @test_nowarn rcs(PECSphere(; radius=2.0), planeWave(; frequency=1e8), [SVector(1.0, 0.0, 0.0)])

    @test_throws ErrorException("The (bistatic) RCS is only defined for a plane-wave excitation (so far).") rcs(
        PECSphere(; radius=1.0), SphericalModeTE(; frequency=f), []
    )
end
