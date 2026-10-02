
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

    @testset "Scattered field" begin

        spHard = HardSphere(; radius=spRadius)
        spSoft = SoftSphere(; radius=spRadius)

        directions = [SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0), SVector(1.0, 0.0, 0.0), normalize(SVector(1.0, 2.0, -0.5))]

        @testset "Boundary conditions" begin

            # the closed-form incident field is combined with the series for the scattered field, so
            # that the boundary conditions verify the coefficients of the incident expansion as well
            for n̂ in directions
                point = spRadius * n̂

                @test abs(field(spHard, ex, PressureNormalGradient([point]))[1]) / abs(field(ex, PressureNormalGradient([point]))[1]) <
                    1e-10
                @test abs(field(spSoft, ex, PressureTrace([point]))[1]) / abs(field(ex, PressureTrace([point]))[1]) < 1e-10
            end
        end

        @testset "Reciprocity" begin

            # the total field is the Green's function of the exterior problem and hence symmetric in
            # the position of the source and of the observer
            for (p₀, p₁) in ((SVector(0.0, 0.0, 3.0), SVector(2.0, 1.0, -1.5)), (SVector(-2.5, 0.5, 0.0), SVector(0.3, -4.0, 1.0)))
                for sp in (spHard, spSoft)
                    @test field(sp, SphericalScattering.Acoustic.monopole(; position=p₀, frequency=f), Pressure([p₁]))[1] ≈
                        field(sp, SphericalScattering.Acoustic.monopole(; position=p₁, frequency=f), Pressure([p₀]))[1] rtol = 1e-12
                end
            end
        end

        @testset "Distant monopole" begin

            # a distant monopole acts as a plane wave travelling in -r̂₀ with amplitude exp(-im κ R₀) / (4π R₀)
            point = SVector(0.3, -0.4, 1.6)

            for (R₀, tol) in ((1e4, 1e-3), (1e6, 1e-5))
                position = R₀ * SVector(0.0, 0.0, 1.0)

                exFar = SphericalScattering.Acoustic.monopole(; position=position, frequency=f)
                exPlane = SphericalScattering.Acoustic.planeWave(; frequency=f, direction=(-normalize(position)))

                @test scatteredfield(spHard, exFar, Pressure([point]))[1] ≈
                    cis(-κ * R₀) / (4π * R₀) * scatteredfield(spHard, exPlane, Pressure([point]))[1] rtol = tol
            end
        end

        @testset "Far field" begin

            r = 1e5 * spRadius

            for sp in (spHard, spSoft), n̂ in directions
                FF = scatteredfield(sp, ex, FarField([n̂]))[1]

                @test scatteredfield(sp, ex, Pressure([r * n̂]))[1] ≈ FF * cis(-κ * r) / r rtol = 1e-4
            end
        end

        @testset "Rotational symmetry" begin

            # the scattered field depends on the geometry alone, not on its orientation in space
            ϑ₀ = 0.6
            r = 2.5

            reference = nothing

            for (axis, perp) in (
                (SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0)),
                (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)),
                (normalize(SVector(1.0, 1.0, 1.0)), normalize(SVector(1.0, -1.0, 0.0))),
            )
                exRot = SphericalScattering.Acoustic.monopole(; position=3.0 * axis, frequency=f)
                p = scatteredfield(spHard, exRot, Pressure([r * (cos(ϑ₀) * axis + sin(ϑ₀) * perp)]))[1]

                isnothing(reference) ? (reference = p) : @test p ≈ reference rtol = 1e-12
            end
        end

        @testset "Provided normals" begin

            # as for the plane wave: n̂ ⋅ ∇p = (n̂ ⋅ r̂) ∂p/∂r + (n̂ ⋅ ϑ̂) 1/a ∂p/∂ϑ, the polar derivative
            # stemming from a central difference of the Dirichlet trace along the sphere
            hϑ = 1e-6
            ê = normalize(r₀)

            function rotated(r̂, δ)
                cosϑ = dot(ê, r̂)
                ê⊥ = normalize(r̂ - cosϑ * ê)
                ϑ = acos(clamp(cosϑ, -1.0, 1.0))

                return cos(ϑ + δ) * ê + sin(ϑ + δ) * ê⊥
            end

            for n̂ in (normalize(SVector(1.0, 2.0, -0.5)), SVector(1.0, 0.0, 0.0))

                cosϑ = dot(ê, n̂)
                ϑ̂ = (cosϑ * n̂ - ê) / sqrt(1 - cosϑ^2)

                gᵣ = scatteredfield(spHard, ex, PressureNormalGradient([spRadius * n̂]))[1]
                dpϑ =
                    (
                        scatteredfield(spHard, ex, PressureTrace([rotated(n̂, hϑ)]))[1] -
                        scatteredfield(spHard, ex, PressureTrace([rotated(n̂, -hϑ)]))[1]
                    ) / (2 * hϑ)

                for w in (0.4, -0.9)
                    ñ = normalize(n̂ + w * ϑ̂)

                    @test scatteredfield(spHard, ex, PressureNormalGradient([spRadius * n̂], [ñ]))[1] ≈
                        dot(ñ, n̂) * gᵣ + dot(ñ, ϑ̂) / spRadius * dpϑ rtol = 1e-7
                end
            end
        end

        @testset "Interior and preconditions" begin

            for sp in (spHard, spSoft)
                @test norm(scatteredfield(sp, ex, Pressure(points_cartNF_inside))) == 0.0
            end

            errInside = ErrorException(
                "The monopole has to be located outside the sphere: its distance from the center is smaller than the radius."
            )
            exInside = SphericalScattering.Acoustic.monopole(; position=SVector(0.0, 0.0, 0.5), frequency=f)

            @test_throws errInside scatteredfield(spHard, exInside, Pressure(points_cartNF))
            @test_throws errInside scatteredfield(spHard, exInside, PressureNormalGradient(points_cartFF))

            errSphere = ErrorException("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")

            @test_throws errSphere scatteredfield(PECSphere(; radius=spRadius), ex, Pressure(points_cartNF))
            @test_throws errSphere scatteredfield(PECSphere(; radius=spRadius), ex, PressureNormalGradient(points_cartFF))
        end
    end

    @testset "Incident and total far field" begin

        directions = [SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0), SVector(1.0, 0.0, 0.0), normalize(SVector(1.0, 2.0, -0.5))]

        @testset "Incident far field" begin

            # --- in contrast to a plane wave, a monopole does possess a far field
            for n̂ in directions
                FF = field(ex, n̂, FarField([n̂]))

                @test FF ≈ cis(κ * dot(n̂, r₀)) / (4π) rtol = 1e-14

                # a point source radiates isotropically: the position enters the phase alone
                @test abs(FF) ≈ 1 / (4π) rtol = 1e-14

                # only the direction of the location is relevant
                @test field(ex, 0.3 * n̂, FarField([n̂])) ≈ FF rtol = 1e-14
                @test field(ex, 1e6 * n̂, FarField([n̂])) ≈ FF rtol = 1e-14
            end

            # --- a monopole in the center radiates the same far field in every direction
            exCenter = SphericalScattering.Acoustic.monopole(; position=SVector(0.0, 0.0, 0.0), frequency=f)

            for n̂ in directions
                @test field(exCenter, n̂, FarField([n̂])) ≈ 1 / (4π) rtol = 1e-14
            end

            # --- the far field is the limit of the pressure with exp(-im κ r) / r removed
            for n̂ in directions, (r, tol) in ((1e4, 1e-2), (1e6, 1e-4))
                @test field(ex, r * n̂, Pressure([r * n̂])) ≈ field(ex, n̂, FarField([n̂])) * cis(-κ * r) / r rtol = tol
            end
        end

        @testset "Total far field" begin

            # the total far field is well defined for a monopole, unlike for a plane wave
            spHard = HardSphere(; radius=spRadius)
            spSoft = SoftSphere(; radius=spRadius)

            r = 1e6 * spRadius

            for sp in (spHard, spSoft), n̂ in directions

                tot = field(sp, ex, FarField([n̂]))[1]

                @test tot ≈ field(ex, n̂, FarField([n̂])) + scatteredfield(sp, ex, FarField([n̂]))[1] rtol = 1e-14

                # and it is the limit of the total pressure
                @test field(sp, ex, Pressure([r * n̂]))[1] ≈ tot * cis(-κ * r) / r rtol = 1e-4
            end
        end

        @testset "Array interface" begin

            FF = field(ex, FarField(points_cartFF))

            @test size(FF) == size(points_cartFF)
            @test all(isfinite, FF)
            @test all(x -> abs(abs(x) - 1 / (4π)) < 1e-15, FF)
            @test FF[1] ≈ field(ex, points_cartFF[1], FarField(points_cartFF)) rtol = 1e-14
        end
    end
end
