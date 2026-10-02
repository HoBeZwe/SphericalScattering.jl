
@testset "Plane wave" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    spHard = HardSphere(; radius=spRadius)
    spSoft = SoftSphere(; radius=spRadius)

    ex = SphericalScattering.Acoustic.planeWave(; frequency=f)

    directions = [SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0), SVector(1.0, 0.0, 0.0), normalize(SVector(1.0, 2.0, -0.5))]

    @testset "Plane wave excitation" begin

        @test ex isa AcousticPlaneWave{Float64,Float64,Float64}
        @test ex.direction ≈ SVector(0.0, 0.0, 1.0)

        # the direction is normalized by the inner constructor
        exSkew = SphericalScattering.Acoustic.planeWave(; frequency=f, direction=SVector(0.0, 3.0, 4.0))
        @test norm(exSkew.direction) ≈ 1.0

        @test_throws ErrorException("missing argument `frequency`") SphericalScattering.Acoustic.planeWave()
    end

    @testset "Incident field" begin

        # the expansion underlying the scattered field has to reproduce the incident plane wave:
        # exp(-im * k * d ⋅ r) = Σ (2n + 1) (-im)ⁿ jₙ(kr) Pₙ(cosϑ)
        point = SVector(1.3, -2.0, 4.1)

        r = norm(point)
        kr = κ * r
        s = sqrt(π / 2 / kr)
        cosϑ = dot(ex.direction, point) / r

        u = ComplexF64(0.0)
        Pn₋₁, Pn = 0.0, 1.0

        for n in 0:60
            u += (2 * n + 1) * (-im)^n * (s * SphericalScattering.besselj(n + 0.5, kr)) * Pn
            Pn₋₁, Pn = Pn, ((2 * n + 1) * cosϑ * Pn - n * Pn₋₁) / (n + 1)
        end

        @test u ≈ field(ex, point, Pressure([point])) rtol = 1e-12
    end

    @testset "Boundary conditions" begin

        # --- sound-hard sphere: the normal velocity, and hence ∂p/∂r, vanishes on the surface
        δ = 1e-9  # offset to stay clear of the surface, where the interior convention takes over
        h = 1e-4  # step size of the one-sided second-order difference quotient

        for n̂ in directions
            p = [field(spHard, ex, Pressure([(spRadius + δ + j * h) * n̂]))[1] for j in 0:2]
            dpdr = (-3 * p[1] + 4 * p[2] - p[3]) / (2 * h)

            @test abs(dpdr) / κ < 1e-5
        end

        # --- sound-soft sphere: the total pressure vanishes on the surface
        for n̂ in directions
            @test abs(field(spSoft, ex, Pressure([(spRadius + δ) * n̂]))[1]) < 1e-7
        end
    end

    @testset "Rayleigh limit" begin

        # for ka ≪ 1 the pressure scattered by a rigid sphere approaches
        # p = A (ka)² (a/r) exp(-im * k * r) (cosϑ/2 - 1/3)
        fl = 1e5
        κl = 2π * fl / c

        exl = SphericalScattering.Acoustic.planeWave(; frequency=fl)

        r = 1e9 * spRadius # the asymptotic form of hₙ⁽²⁾ requires kr ≫ 1 as well

        for cosϑ in (1.0, 0.5, 0.0, -0.5, -1.0)
            point = r * SVector(sqrt(1 - cosϑ^2), 0.0, cosϑ)

            p = scatteredfield(spHard, exl, Pressure([point]))[1]
            pRayleigh = (κl * spRadius)^2 * (spRadius / r) * cis(-κl * r) * (cosϑ / 2 - 1 / 3)

            @test p ≈ pRayleigh rtol = 1e-3
        end
    end

    @testset "Rotational symmetry" begin

        # the scattered pressure depends on the observation point solely via r and d ⋅ r̂
        ϑ₀ = 0.7
        r = 3.3

        ex₁ = SphericalScattering.Acoustic.planeWave(; frequency=f, direction=SVector(0.0, 0.0, 1.0))
        p₁ = scatteredfield(spHard, ex₁, Pressure([r * SVector(sin(ϑ₀), 0.0, cos(ϑ₀))]))[1]

        ex₂ = SphericalScattering.Acoustic.planeWave(; frequency=f, direction=SVector(1.0, 0.0, 0.0))
        p₂ = scatteredfield(spHard, ex₂, Pressure([r * SVector(cos(ϑ₀), sin(ϑ₀), 0.0)]))[1]

        d₃ = normalize(SVector(1.0, 1.0, 1.0))
        e₃ = normalize(cross(d₃, SVector(0.0, 0.0, 1.0)))
        ex₃ = SphericalScattering.Acoustic.planeWave(; frequency=f, direction=d₃)
        p₃ = scatteredfield(spHard, ex₃, Pressure([r * (cos(ϑ₀) * d₃ + sin(ϑ₀) * e₃)]))[1]

        @test p₁ ≈ p₂ rtol = 1e-12
        @test p₁ ≈ p₃ rtol = 1e-12
    end

    @testset "Interior" begin

        for sp in (spHard, spSoft)
            @test norm(scatteredfield(sp, ex, Pressure(points_cartNF_inside))) == 0.0
            @test norm(field(sp, ex, Pressure(points_cartNF_inside))) == 0.0
        end
    end

    @testset "Array interface" begin

        for sp in (spHard, spSoft)
            F = scatteredfield(sp, ex, Pressure(points_cartNF))

            @test size(F) == size(points_cartNF)
            @test all(isfinite, F)

            # the array and the single-point methods have to agree
            @test F[1] ≈ scatteredfield(sp, ex, points_cartNF[1], Pressure(points_cartNF))
        end
    end

    @testset "Far field" begin

        # --- the far field is defined as lim(r → ∞) r exp(im * κ * r) p, so that the pressure has to
        #     approach FF exp(-im * κ * r) / r with an error of the order 1 / (κr)
        r = 1e6 * spRadius

        for sp in (spHard, spSoft), n̂ in directions
            FF = scatteredfield(sp, ex, FarField([n̂]))[1]
            p  = scatteredfield(sp, ex, Pressure([r * n̂]))[1]

            @test p ≈ FF * cis(-κ * r) / r rtol = 1e-5
        end

        # --- the far field is determined by the direction of observation alone: it does not depend on
        #     the distance and, in particular, does not vanish inside the sphere
        for n̂ in directions
            FF = scatteredfield(spHard, ex, FarField([n̂]))[1]

            @test !iszero(FF)
            @test scatteredfield(spHard, ex, FarField([0.5 * spRadius * n̂]))[1] ≈ FF rtol = 1e-12
            @test scatteredfield(spHard, ex, FarField([1e6 * spRadius * n̂]))[1] ≈ FF rtol = 1e-12
        end

        # --- for ka ≪ 1 the far field of a rigid sphere approaches A k² a³ (cosϑ/2 - 1/3)
        fl = 1e5
        κl = 2π * fl / c

        exl = SphericalScattering.Acoustic.planeWave(; frequency=fl)

        for cosϑ in (1.0, 0.5, 0.0, -0.5, -1.0)
            point = SVector(sqrt(1 - cosϑ^2), 0.0, cosϑ)

            FF = scatteredfield(spHard, exl, FarField([point]))[1]

            @test FF ≈ κl^2 * spRadius^3 * (cosϑ / 2 - 1 / 3) rtol = 1e-4
        end

        # --- array interface: the far-field points lie on the surface of the sphere, as for the
        #     electromagnetic excitations
        for sp in (spHard, spSoft)
            FF = scatteredfield(sp, ex, FarField(points_cartFF))

            @test size(FF) == size(points_cartFF)
            @test all(isfinite, FF)
            @test FF[1] ≈ scatteredfield(sp, ex, points_cartFF[1], FarField(points_cartFF))
        end
    end

    @testset "Traces" begin

        # --- spherical Bessel and Hankel functions, for the closed forms below
        jn(n, x)  = sqrt(π / 2 / x) * SphericalScattering.besselj(n + 0.5, x)
        hn(n, x)  = sqrt(π / 2 / x) * SphericalScattering.hankelh2(n + 0.5, x)
        dhn(n, x) = hn(n - 1, x) - (n + 1) / x * hn(n, x)

        ka = κ * spRadius

        @testset "Incident traces" begin

            for n̂ in directions
                point = spRadius * n̂

                # the Dirichlet trace is the incident pressure itself
                @test field(ex, PressureTrace([point]))[1] == field(ex, Pressure([point]))[1]

                # the Neumann trace is -im * κ * (d ⋅ n̂) times the incident pressure
                cosϑ = dot(ex.direction, n̂)
                @test field(ex, PressureNormalGradient([point]))[1] ≈ -im * κ * cosϑ * field(ex, Pressure([point]))[1] rtol = 1e-12
            end
        end

        @testset "Scattered traces" begin

            h = 1e-5

            for sp in (spHard, spSoft), n̂ in directions

                # the Dirichlet trace is the scattered pressure evaluated at r = a
                @test scatteredfield(sp, ex, PressureTrace([n̂]))[1] == scatteredfield(sp, ex, Pressure([spRadius * n̂]))[1]

                # the Neumann trace is the radial derivative of the scattered pressure at r = a
                g = scatteredfield(sp, ex, PressureNormalGradient([n̂]))[1]
                p = [scatteredfield(sp, ex, Pressure([(spRadius + 1e-9 + j * h) * n̂]))[1] for j in 0:2]

                @test abs(g - (-3 * p[1] + 4 * p[2] - p[3]) / (2 * h)) / κ < 1e-6
            end
        end

        @testset "Boundary conditions" begin

            # the total traces fulfill the boundary conditions exactly
            for n̂ in directions
                @test abs(field(spHard, ex, PressureNormalGradient([n̂]))[1]) / κ < 1e-10
                @test abs(field(spSoft, ex, PressureTrace([n̂]))[1]) < 1e-10
            end
        end

        @testset "Closed forms of the total traces" begin

            # eliminating bₙ with the Wronskian jₙ(x) yₙ′(x) - jₙ′(x) yₙ(x) = 1/x² yields
            # p = -im / (ka)² Σ (2n+1) (-im)ⁿ Pₙ / hₙ′(ka)      on a sound-hard sphere and
            # ∂p/∂n = im κ / (ka)² Σ (2n+1) (-im)ⁿ Pₙ / hₙ(ka)  on a sound-soft sphere
            for n̂ in directions

                cosϑ = dot(ex.direction, n̂)

                uHard, uSoft = ComplexF64(0.0), ComplexF64(0.0)
                Pn₋₁, Pn = 0.0, 1.0

                for n in 0:60
                    uHard += (2 * n + 1) * (-im)^n * Pn / dhn(n, ka)
                    uSoft += (2 * n + 1) * (-im)^n * Pn / hn(n, ka)
                    Pn₋₁, Pn = Pn, ((2 * n + 1) * cosϑ * Pn - n * Pn₋₁) / (n + 1)
                end

                @test field(spHard, ex, PressureTrace([n̂]))[1] ≈ -im / ka^2 * uHard rtol = 1e-10
                @test field(spSoft, ex, PressureNormalGradient([n̂]))[1] ≈ im * κ / ka^2 * uSoft rtol = 1e-10
            end
        end

        @testset "Evaluation on the surface" begin

            # the scattered traces are evaluated on the sphere, so that locations of a faceted surface
            # mesh, which do not lie exactly on the sphere, yield the same values
            for n̂ in directions, quantity in (PressureTrace, PressureNormalGradient)

                ref = scatteredfield(spHard, ex, quantity([spRadius * n̂]))[1]

                # the normals are re-derived from the scaled locations, so rounding at the last
                # digit is admissible
                for scale in (0.93, 0.999, 1.05)
                    @test scatteredfield(spHard, ex, quantity([scale * spRadius * n̂]))[1] ≈ ref rtol = 1e-14
                end
            end
        end

        @testset "Provided normals" begin

            # --- the normals default to the normalized locations
            q = PressureNormalGradient(points_cartFF)

            @test q.normals ≈ map(normalize, points_cartFF)
            @test size(q.normals) == size(points_cartFF)

            # --- providing the radial normals explicitly reproduces the default
            @test field(ex, PressureNormalGradient(points_cartFF, map(normalize, points_cartFF))) ≈
                field(ex, PressureNormalGradient(points_cartFF)) rtol = 1e-14

            # --- the normals are normalized, hence their length is irrelevant
            @test field(ex, PressureNormalGradient(points_cartFF, map(p -> 7.3 * p, points_cartFF))) ≈
                field(ex, PressureNormalGradient(points_cartFF)) rtol = 1e-14

            # --- a normal that is not radial: compare against the directional derivative of the
            #     incident pressure, as is relevant for the facets of a surface mesh
            h = 1e-6

            for n̂ in directions
                t̂ = normalize(cross(n̂, SVector(0.0, 0.0, 1.0)) + cross(n̂, SVector(1.0, 0.0, 0.0)))
                ñ = normalize(n̂ + 0.3 * t̂)   # tilted away from r̂, as for a flat facet

                point = spRadius * n̂

                g = field(ex, PressureNormalGradient([point], [ñ]))[1]
                fd = (field(ex, Pressure([point + h * ñ]))[1] - field(ex, Pressure([point - h * ñ]))[1]) / (2 * h)

                @test abs(g - fd) / κ < 1e-7

                # the analytic expression, for good measure
                @test g ≈ -im * κ * dot(ex.direction, ñ) * field(ex, Pressure([point]))[1] rtol = 1e-12
            end

            # --- the number of normals has to match the number of locations
            @test_throws ErrorException("The number of provided normal vectors does not match the number of locations.") PressureNormalGradient(
                points_cartFF, [SVector(0.0, 0.0, 1.0)]
            )

        end

        @testset "Provided normals for the scattered field" begin

            # n̂ ⋅ ∇p = (n̂ ⋅ r̂) ∂p/∂r + (n̂ ⋅ ϑ̂) 1/a ∂p/∂ϑ, where ∂p/∂ϑ is obtained from a central
            # difference of the Dirichlet trace along the sphere, that is, at a constant radius
            hϑ = 1e-6

            # rotate a surface point within the (d, r̂) plane, staying exactly on the sphere
            function rotated(r̂, δ)
                cosϑ = dot(ex.direction, r̂)
                ê = normalize(r̂ - cosϑ * ex.direction)
                ϑ = acos(clamp(cosϑ, -1.0, 1.0))

                return cos(ϑ + δ) * ex.direction + sin(ϑ + δ) * ê
            end

            polar(r̂) = (dot(ex.direction, r̂) * r̂ - ex.direction) / sqrt(1 - dot(ex.direction, r̂)^2)

            # ϑ̂ is defined away from the poles only
            offAxis = [normalize(SVector(1.0, 2.0, -0.5)), SVector(1.0, 0.0, 0.0), normalize(SVector(-0.3, 0.4, 0.6))]

            for sp in (spHard, spSoft), n̂ in offAxis

                ϑ̂ = polar(n̂)

                gᵣ = scatteredfield(sp, ex, PressureNormalGradient([spRadius * n̂]))[1]
                dpϑ =
                    (
                        scatteredfield(sp, ex, PressureTrace([rotated(n̂, hϑ)]))[1] -
                        scatteredfield(sp, ex, PressureTrace([rotated(n̂, -hϑ)]))[1]
                    ) / (2 * hϑ)

                for w in (0.3, 1.0, -0.7)
                    ñ = normalize(n̂ + w * ϑ̂)

                    @test scatteredfield(sp, ex, PressureNormalGradient([spRadius * n̂], [ñ]))[1] ≈
                        dot(ñ, n̂) * gᵣ + dot(ñ, ϑ̂) / spRadius * dpϑ rtol = 1e-7
                end

                # a purely tangential normal picks out the polar derivative alone
                @test scatteredfield(sp, ex, PressureNormalGradient([spRadius * n̂], [ϑ̂]))[1] ≈ dpϑ / spRadius rtol = 1e-7
            end

            # --- on a rigid sphere ∂p_tot/∂r vanishes, so an arbitrary normal picks out the
            #     tangential part of the total gradient alone
            for n̂ in offAxis

                ϑ̂ = polar(n̂)
                ñ = normalize(n̂ + 0.5 * ϑ̂)

                dpϑ =
                    (field(spHard, ex, PressureTrace([rotated(n̂, hϑ)]))[1] - field(spHard, ex, PressureTrace([rotated(n̂, -hϑ)]))[1]) /
                    (2 * hϑ)

                @test field(spHard, ex, PressureNormalGradient([spRadius * n̂], [ñ]))[1] ≈ dot(ñ, ϑ̂) / spRadius * dpϑ rtol = 1e-7
            end

            # --- the total pressure vanishes on a sound-soft sphere, hence so does its polar
            #     derivative: the total gradient is radial for every normal
            for n̂ in offAxis

                ϑ̂ = polar(n̂)
                gᵣ = field(spSoft, ex, PressureNormalGradient([spRadius * n̂]))[1]

                for w in (0.4, -1.3)
                    ñ = normalize(n̂ + w * ϑ̂)

                    @test field(spSoft, ex, PressureNormalGradient([spRadius * n̂], [ñ]))[1] ≈ dot(ñ, n̂) * gᵣ rtol = 1e-11
                end
            end
        end

        @testset "Rayleigh limit" begin

            # for ka ≪ 1 the surface pressure of a rigid sphere approaches A (1 - 3/2 im ka cosϑ)
            fl = 1e5
            κl = 2π * fl / c

            exl = SphericalScattering.Acoustic.planeWave(; frequency=fl)

            for cosϑ in (1.0, 0.5, 0.0, -1.0)
                point = SVector(sqrt(1 - cosϑ^2), 0.0, cosϑ)

                @test field(spHard, exl, PressureTrace([point]))[1] ≈ 1 - 1.5im * κl * spRadius * cosϑ rtol = 1e-4
            end
        end

        @testset "Array interface" begin

            for sp in (spHard, spSoft), quantity in (PressureTrace, PressureNormalGradient)
                F = scatteredfield(sp, ex, quantity(points_cartFF))
                G = field(sp, ex, quantity(points_cartFF))

                @test size(F) == size(points_cartFF)
                @test size(G) == size(points_cartFF)
                @test all(isfinite, F)
                @test all(isfinite, G)
            end
        end

        @testset "Unsupported spheres" begin

            err = ErrorException("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")

            @test_throws err scatteredfield(PECSphere(; radius=spRadius), ex, PressureTrace(points_cartFF))
            @test_throws err scatteredfield(PECSphere(; radius=spRadius), ex, PressureNormalGradient(points_cartFF))
        end
    end

    @testset "Undefined far fields" begin

        @test_throws ErrorException("The far-field of a plane wave is not defined.") field(ex, FarField(points_cartFF))

        @test_throws ErrorException("The total far-field for a plane-wave excitation is not defined") field(
            spHard, ex, FarField(points_cartFF)
        )
    end

    @testset "Unsupported spheres" begin

        err = ErrorException("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")

        @test_throws err scatteredfield(PECSphere(; radius=spRadius), ex, Pressure(points_cartNF))
        @test_throws err scatteredfield(PECSphere(; radius=spRadius), ex, FarField(points_cartFF))
    end
end
