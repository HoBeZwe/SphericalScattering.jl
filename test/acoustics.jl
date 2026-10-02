
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
