
@testitem "Spherical coordinates and rotations" setup = [Setup] begin

    # ----- sph2cart
    vec = SVector(1.0, π / 4, π / 4)
    xyz = SphericalScattering.sph2cart(vec)

    ẑ = SVector(0.0, 0.0, 1.0)
    rθϕ = SphericalScattering.cart2sph(-ẑ)

    @test rθϕ[2] ≈ π

    @test xyz[1] ≈ 0.5
    @test xyz[2] ≈ 0.5
    @test xyz[3] ≈ 1 / √2

    # ----- convertCartesian2Spherical
    F_cart    = SVector(1.0, 1.0, 1.0)
    point_sph = SVector(1, π / 2, π / 2)
    F_sph     = SphericalScattering.convertCartesian2Spherical(F_cart, point_sph)

    @test F_sph[1] ≈ +1.0
    @test F_sph[2] ≈ -1.0
    @test F_sph[3] ≈ -1.0

    # ----- rotation matrix corner case
    orient = normalize(SVector(0.0, 0.0, -1.0))
    ex = HertzianDipole(; frequency=1e8, orientation=orient, position=2 * orient)

    @test_nowarn SphericalScattering.rotationMatrix(ex)
end


@testitem "Oblate spheroidal coordinates" setup = [Setup] begin

    f = 1.3   # semifocal distance

    @testset "Coordinate transforms" begin

        points = [
            SVector(1.0, 2.0, -0.5),
            SVector(0.0, 0.0, 4.0),
            SVector(3.0, -1.0, 0.0),
            SVector(-0.2, 0.1, 0.05),
            SVector(0.0, 0.0, 0.0),
            SVector(f, 0.0, 0.0),      # the rim of the disc
        ]

        # --- the transforms are inverse to each other
        for point in points
            @test SphericalScattering.obl2cart(SphericalScattering.cart2obl(point, f), f) ≈ point atol = 1e-14
        end

        # --- the coordinates lie in their intended ranges
        for point in points
            ξ, η, φ = SphericalScattering.cart2obl(point, f)

            @test ξ >= 0
            @test -1 <= η <= 1
            @test -π <= φ <= π
        end

        # --- ξ = const is an oblate spheroid with the expected semi-axes
        for ξ in (0.0, 0.3, 1.0, 5.0)
            a = f * sqrt(1 + ξ^2)
            b = f * ξ

            for η in (-1.0, -0.4, 0.0, 0.7), φ in (0.0, 1.2, 4.0)
                x, y, z = SphericalScattering.obl2cart(SVector(ξ, η, φ), f)

                if iszero(ξ)
                    @test abs(z) < 1e-15          # the disc lies in the plane z = 0
                    @test hypot(x, y) <= f + 1e-15 # and has radius f
                else
                    @test (x^2 + y^2) / a^2 + z^2 / b^2 ≈ 1 rtol = 1e-14
                end
            end
        end

        # --- the rim and the center of the disc
        @test SphericalScattering.cart2obl(SVector(f, 0.0, 0.0), f)[1] ≈ 0.0 atol = 1e-14
        @test SphericalScattering.cart2obl(SVector(f, 0.0, 0.0), f)[2] ≈ 0.0 atol = 1e-7
        @test SphericalScattering.cart2obl(SVector(0.0, 0.0, 0.0), f)[1] ≈ 0.0 atol = 1e-14
        @test abs(SphericalScattering.cart2obl(SVector(0.0, 0.0, 0.0), f)[2]) ≈ 1.0 atol = 1e-14
    end

    @testset "Metric coefficients" begin

        # the metric coefficients have to be the lengths of the coordinate tangent vectors, which are
        # obtained independently from central differences of the forward transform
        h = 1e-6

        for ξ in (0.2, 1.0, 3.0), η in (-0.7, -0.1, 0.4), φ in (0.0, 2.0)

            point = SVector(ξ, η, φ)
            metric = SphericalScattering.oblateMetric(point, f)

            for (ind, δ) in enumerate((SVector(h, 0.0, 0.0), SVector(0.0, h, 0.0), SVector(0.0, 0.0, h)))
                tangent = (SphericalScattering.obl2cart(point + δ, f) - SphericalScattering.obl2cart(point - δ, f)) / (2 * h)

                @test norm(tangent) ≈ metric[ind] rtol = 1e-8
            end
        end

        # on the disc the normal derivative is singular at the rim, since h_ξ = f |η|
        @test SphericalScattering.oblateMetric(SVector(0.0, 0.5, 0.0), f)[1] ≈ f * 0.5 rtol = 1e-14
        @test SphericalScattering.oblateMetric(SVector(0.0, 0.0, 0.0), f)[1] ≈ 0.0 atol = 1e-15
    end
end


@testitem "Prolate spheroidal coordinates" setup = [Setup] begin

    SS = SphericalScattering

    f = 0.8

    # --- the transforms are inverse to each other, the coordinates lie in their ranges
    for x in (-1.3, 0.0, 0.4, 2.0), y in (0.0, -0.7, 0.9), z in (-1.1, 0.0, 0.3, 0.79, 1.7)
        point = SVector(x, y, z)
        ξ, η, φ = SS.cart2prol(point, f)

        @test SS.prol2cart(SVector(ξ, η, φ), f) ≈ point atol = 1e-14
        @test ξ >= 1
        @test -1 <= η <= 1
    end

    # --- on the axis beyond the foci η = ±1 holds exactly, on the segment between them ξ = 1 and η = z / f
    @test SS.cart2prol(SVector(0.0, 0.0, 1.3), f)[2] == 1.0
    @test SS.cart2prol(SVector(0.0, 0.0, -1.3), f)[2] == -1.0
    @test SS.cart2prol(SVector(0.0, 0.0, 0.5), f)[1:2] ≈ SVector(1.0, 0.5 / f) rtol = 1e-15

    # --- in the plane z = 0 the sign of a vanishing z does not matter
    @test SS.cart2prol(SVector(0.3, -0.2, -0.0), f) == SS.cart2prol(SVector(0.3, -0.2, 0.0), f)

    # --- ξ = const is a prolate spheroid with the semi-axes f √(ξ² - 1) and f ξ
    for ξ in (1.05, 1.5, 3.0), η in (-0.9, -0.2, 0.0, 0.6, 1.0), φ in (0.0, 1.3)
        x, y, z = SS.prol2cart(SVector(ξ, η, φ), f)
        @test (x^2 + y^2) / (f^2 * (ξ^2 - 1)) + z^2 / (f * ξ)^2 ≈ 1 rtol = 1e-14
    end

    # --- the metric coefficients are the lengths of the tangent vectors, the basis is orthonormal
    h = 1e-6
    for (ξ, η, φ) in ((1.3, 0.4, 0.7), (2.1, -0.8, 2.5), (1.02, 0.1, -1.2))
        q = SVector(ξ, η, φ)
        metric = SS.prolateMetric(q, f)
        basis = SS.prolateBasis(q, f)

        for (ind, e) in enumerate((SVector(h, 0, 0), SVector(0, h, 0), SVector(0, 0, h)))
            tangent = (SS.prol2cart(q + e, f) - SS.prol2cart(q - e, f)) / (2h)

            @test norm(tangent) ≈ metric[ind] rtol = 1e-8
            @test tangent / norm(tangent) ≈ basis[ind] atol = 1e-8
        end

        @test abs(dot(basis[1], basis[2])) < 1e-14
        @test abs(dot(basis[1], basis[3])) < 1e-14
        @test abs(dot(basis[2], basis[3])) < 1e-14
    end
end
