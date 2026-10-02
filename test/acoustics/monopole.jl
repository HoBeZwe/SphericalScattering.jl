
@testset "Monopole" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    r₀ = SVector(0.0, 0.0, 3.0)   # well clear of the evaluation points, the field being singular there

    ex = SphericalScattering.Acoustic.monopole(; position=r₀, frequency=f)

    points = [SVector(1.0, 2.0, -0.5), SVector(0.0, 0.0, 1.0), SVector(-1.3, 0.4, 0.7), SVector(2.0, -1.0, 5.0)]
    normals = [normalize(SVector(1.0, -0.4, 0.3)), SVector(0.0, 1.0, 0.0), normalize(SVector(-1.0, 1.0, 2.0)), SVector(1.0, 0.0, 0.0)]

    @testset "Monopole excitation" begin

        @test ex isa AcousticMonopole{Float64,Float64,Float64}
        @test ex.position == r₀

        @test_throws ErrorException("missing argument `position`") SphericalScattering.Acoustic.monopole(; frequency=f)
        @test_throws ErrorException("missing argument `frequency`") SphericalScattering.Acoustic.monopole(; position=r₀)
    end

    @testset "Incident field" begin

        # --- the pressure is the free-space Green's function scaled by the amplitude
        for point in points
            R = norm(point - r₀)

            @test field(ex, point, Pressure([point])) ≈ cis(-κ * R) / (4π * R) rtol = 1e-14
        end

        # --- the amplitude enters linearly
        exScaled = SphericalScattering.Acoustic.monopole(; position=r₀, frequency=f, amplitude=2.5)

        @test field(exScaled, Pressure(points)) ≈ 2.5 * field(ex, Pressure(points)) rtol = 1e-14

        # --- the pressure has to fulfill the Helmholtz equation ∇²p + κ²p = 0, which a point source
        #     without the factor 1/R would not
        h = 1e-5

        for point in points
            p₀ = field(ex, point, Pressure([point]))

            laplacian = ComplexF64(0.0)
            for i in 1:3
                eᵢ = SVector(ntuple(j -> j == i ? 1.0 : 0.0, 3))
                laplacian +=
                    (field(ex, point + h * eᵢ, Pressure([point])) - 2 * p₀ + field(ex, point - h * eᵢ, Pressure([point]))) / h^2
            end

            @test abs(laplacian + κ^2 * p₀) / (κ^2 * abs(p₀)) < 1e-4
        end

        # --- the incident field is suppressed inside a sphere for the total field
        @test iszero(field(ex, Pressure([SVector(0.0, 0.0, 0.5)]); zeroRadius=1.0)[1])
    end

    @testset "Traces" begin

        # --- the Dirichlet trace is the pressure itself
        @test field(ex, PressureTrace(points)) == field(ex, Pressure(points))

        # --- the normals default to the normalized locations
        @test PressureNormalGradient(points).normals == normalize.(points)
        @test field(ex, PressureNormalGradient(points))[1] ==
            field(ex, points[1], normalize(points[1]), PressureNormalGradient(points))

        # --- the Neumann trace against the analytic expression and a directional finite difference
        h = 1e-6

        for (point, n̂) in zip(points, normals)
            R = norm(point - r₀)
            R̂ = (point - r₀) / R

            g = field(ex, point, n̂, PressureNormalGradient([point], [n̂]))

            @test g ≈ -(1 + im * κ * R) * cis(-κ * R) / (4π * R^2) * dot(n̂, R̂) rtol = 1e-14

            fd = (field(ex, point + h * n̂, Pressure([point])) - field(ex, point - h * n̂, Pressure([point]))) / (2 * h)

            @test g ≈ fd rtol = 1e-7
        end

        # --- the array method employs the normals of the quantity
        @test field(ex, PressureNormalGradient(points, normals)) ≈
            [field(ex, points[i], normals[i], PressureNormalGradient(points, normals)) for i in eachindex(points)] rtol = 1e-14
    end

    @testset "Limiting cases" begin

        # --- for κR ≫ 1 the Neumann trace approaches the local plane-wave result -im κ (R̂ ⋅ n̂) p
        n̂ = normalize(SVector(0.3, 1.0, -0.2))

        for (R, tol) in ((1e3, 1e-3), (1e5, 1e-5))
            point = r₀ + R * normalize(SVector(1.0, 1.0, 1.0))

            g = field(ex, point, n̂, PressureNormalGradient([point], [n̂]))

            @test g ≈ -im * κ * dot(normalize(point - r₀), n̂) * field(ex, point, Pressure([point])) rtol = 10 * tol
        end

        # --- for κR ≪ 1 the static limit of a point source remains
        exStatic = SphericalScattering.Acoustic.monopole(; position=r₀, frequency=1e1)

        point = points[1]
        R = norm(point - r₀)
        R̂ = (point - r₀) / R

        @test field(exStatic, point, Pressure([point])) ≈ 1 / (4π * R) rtol = 1e-5
        @test field(exStatic, point, normals[1], PressureNormalGradient([point], [normals[1]])) ≈ -dot(normals[1], R̂) / (4π * R^2) rtol =
            1e-10
    end
end
