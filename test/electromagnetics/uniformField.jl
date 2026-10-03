
@testsnippet UniformFieldSetup begin

    # μ0 and ε0 are exported by the package, hence not redefined here
    dir = SVector(1.0, 5.0, -3.0)
    ex = UniformField(; direction=dir, amplitude=norm(dir))
end

@testitem "Incident field" setup = [Setup, UniformFieldSetup] begin
    # define an observation point
    point_cart = [SVector(2.0, 2.0, 3.2)]

    E = field(ex, ElectricField(point_cart))
    @test E[1][1] ≈ 1.0
    @test E[1][2] ≈ 5.0
    @test E[1][3] ≈ -3.0
end

@testitem "Scattered fields" setup = [Setup, UniformFieldSetup] begin
    @testset "Dielectric sphere" begin
        # define scatterer: dielectric sphere
        sp = DielectricSphere(; radius=1.0, filling=Medium(ε0 * 5, μ0))

        # define observation points inside and outside of the sphere
        point_cart = [SVector(0.0, 0.0, 0.5), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

        @test Φ[1] ≈ -6 / 7
        @test Φ[2] ≈ 1 / 7

        @test E[1][1] ≈ -4 / 7
        @test E[1][2] ≈ -20 / 7
        @test E[1][3] ≈ 12 / 7

        @test E[2][1] ≈ 1 / 7
        @test E[2][2] ≈ -5 / 14
        @test E[2][3] ≈ 3 / 14
    end
    @testset "PEC sphere" begin
        # define scatterer: PEC sphere
        sp = PECSphere(1.0)

        # define observation points inside and outside of the sphere
        point_cart = [SVector(0.0, 0.0, 0.5), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

        @test Φ[1] ≈ -3 / 2
        @test Φ[2] ≈ 1 / 4

        @test E[1][1] ≈ -1
        @test E[1][2] ≈ -5
        @test E[1][3] ≈ 3

        @test E[2][1] ≈ 1 / 4
        @test E[2][2] ≈ -5 / 8
        @test E[2][3] ≈ 3 / 8
    end
    @testset "Layered sphere" begin
        # define scatterer: layered dielectric sphere
        sp = LayeredSphere(; radii=SVector(0.5, 1.0), filling=SVector(Medium(5ε0, μ0), Medium(3ε0, μ0)))

        # define observation points in both layers and outside of the sphere
        point_cart = [SVector(0.25, 0.0, 0.0), SVector(0.75, 0.0, 0.0), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

        # Sihvola&Lindell 1988, (10) - (12)
        @test Φ[1] ≈ 29 / 224
        @test Φ[2] ≈ 223 / 672
        @test Φ[3] ≈ 95 / 896

        @test E[1][1] ≈ -29 / 56
        @test E[1][2] ≈ -145 / 56
        @test E[1][3] ≈ 87 / 56

        @test E[2][1] ≈ -25 / 72
        @test E[2][2] ≈ -1115 / 504
        @test E[2][3] ≈ 223 / 168

        @test E[3][1] ≈ 95 / 896
        @test E[3][2] ≈ -475 / 1792
        @test E[3][3] ≈ 285 / 1792
    end
    @testset "Layered sphere PEC" begin
        # define scatterer: layered sphere PEC
        sp = LayeredSpherePEC(; radii=SVector(0.5, 1.0), filling=SVector(Medium(5ε0, μ0)))

        # define observation points in both layers and outside of the sphere
        point_cart = [SVector(0.25, 0.0, 0.0), SVector(0.75, 0.0, 0.0), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

        @test Φ[1] ≈ 1 / 4
        @test Φ[2] ≈ 53 / 96
        @test Φ[3] ≈ 43 / 256

        @test E[1][1] ≈ -1
        @test E[1][2] ≈ -5
        @test E[1][3] ≈ 3

        @test E[2][1] ≈ -29 / 72
        @test E[2][2] ≈ -265 / 72
        @test E[2][3] ≈ 53 / 24

        @test E[3][1] ≈ 43 / 256
        @test E[3][2] ≈ -215 / 512
        @test E[3][3] ≈ 129 / 512


        # more than one dielectric layer
        sp = LayeredSpherePEC(; radii=SVector(0.5, 0.75, 1.0), filling=SVector(Medium(3ε0, μ0), Medium(2ε0, μ0)))

        # define observation points in both layers and outside of the sphere
        point_cart = [SVector(0.25, 0.0, 0.0), SVector(0.6, 0.0, 0.0), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

        Φ[1] ≈ 1 / 4
        Φ[2] ≈ 13287 / 28697
        Φ[3] ≈ 25405 / 229576

        E[1][1] ≈ -1
        E[1][2] ≈ -5
        E[1][3] ≈ 3
    end
    @testset "Dielectric sphere with thin impedance layer" begin
        # Free-space permittivity
        ε₀ = 8.854187812813e-18 # F µm⁻¹ - > 8.854.. * e-12 F m⁻¹
        E₀ = 1e-6
        # Radius of cell
        R_c = 1 # µm
        Δ = 5e-3 #5e-3 # µm Membrane thickness

        # Permittivity and conductivity outside (suspension)
        ε_s = 80 #
        σ_s = 1e-6 # S / µm --> 1 S / m

        #
        ε_m = 10 # Permittivity
        σ_m = 1e-12 # S / µm --> 1e-6 S / m

        σ_c = 1e-6

        f = 1e7
        ω = 2 * π * f

        # Assume exp(im*ω*t)-time dependency
        md_s = Medium(ε_s * ε₀ + im * (σ_s / ω), σ_s)
        md_m = Medium(ε_m * ε₀ + im * (σ_m / ω), σ_m)
        md_c = Medium(ε_s * ε₀ + im * (σ_c / ω), σ_m)

        #
        sp = LayeredSphere(;
            radii=SVector(R_c - Δ, R_c),
            filling=SVector(md_c, md_m), # From inner to outer layer
        )


        spj = DielectricSphereThinImpedanceLayer(; radius=R_c, thickness=Δ, thinlayer=md_m, filling=md_c)

        # We compare against Jones, 1995, Appendix, he uses E₀ẑ
        # So potential points (i.e., grows) into ẑ direction
        potential_direction = dir

        ex = UniformField(; amplitude=E₀, direction=(-potential_direction), embedding=md_s)

        Φsca_ana_3l(pts) = scatteredfield(sp, ex, ScalarPotential(pts))
        Φtot_ana_3l(pts) = field(sp, ex, ScalarPotential(pts))
        ∇Φsca_ana_3l(pts) = scatteredfield(sp, ex, ElectricField(pts)) # I want the gradient, not the electric field

        Φsca_ana_app(pts) = scatteredfield(spj, ex, ScalarPotential(pts))
        Φtot_ana_app(pts) = field(spj, ex, ScalarPotential(pts))
        ∇Φsca_ana_app(pts) = scatteredfield(spj, ex, ElectricField(pts)) # I want the gradient, not the electric field

        @test norm(Φsca_ana_app(points_cartNF) - Φsca_ana_3l(points_cartNF)) / norm(Φsca_ana_3l(points_cartNF)) < 0.007
        @test norm(Φsca_ana_app(points_cartNF_inside) - Φsca_ana_3l(points_cartNF_inside)) / norm(Φsca_ana_3l(points_cartNF_inside)) <
            0.02
        @test norm(∇Φsca_ana_app(points_cartNF) - ∇Φsca_ana_3l(points_cartNF)) / norm(∇Φsca_ana_3l(points_cartNF)) < 0.007

        # Jump of potential should be comparable to the exact model
        ΔΦ_3l(pts) = Φtot_ana_3l(pts ./ norm.(pts) .* ((R_c - Δ))) .- Φtot_ana_3l(pts ./ norm.(pts) .* R_c)
        ΔΦ(pts) = scatteredfield(spj, ex, ScalarPotentialJump(pts))

        tmp_approx_relerror(pts) = norm(ΔΦ(pts) - ΔΦ_3l(pts)) / norm(ΔΦ_3l(pts))

        tmp_approx_relerror(points_cartNF) < 0.0015

        ε∇Φsca_ana_app(pts) = scatteredfield(spj, ex, DisplacementField(pts))

        Esca(pts) = scatteredfield(spj, ex, ElectricField(pts))
        𝒏 = points_cartFF ./ norm.(points_cartFF)
        absdiff = dot.(𝒏, ex.embedding.ε * Esca(points_cartFF .* 1.01) - ε∇Φsca_ana_app(points_cartFF .* 0.99))

        # Check that normal component of D-field is continuous
        @test norm(absdiff) / norm(dot.(𝒏, ε∇Φsca_ana_app(points_cartFF .* 0.99))) < 0.03


        @test norm(ε∇Φsca_ana_app(points_cartNF * 2.0) - ex.embedding.ε * Esca(points_cartNF * 2.0)) /
              abs(ex.embedding.ε * norm(Esca(points_cartNF * 2.0))) < 0.007
    end
end
