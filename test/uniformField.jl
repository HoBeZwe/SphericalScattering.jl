using StaticArrays

const μ0 = 4pi * 1e-7        # default permeability
const ε0 = 8.8541878176e-12  # default permittivity
ex = UniformField(; direction=SVector(1.0, 5.0, -3.0))

@testset "Incident field" begin
    # define an observation point
    point_cart = [SVector(2.0, 2.0, 3.2)]

    E = field(ex, ElectricField(point_cart))
    @test E[1][1] == 1.0
    @test E[1][2] == 5.0
    @test E[1][3] == -3.0
end

@testset "Scattered fields" begin
    @testset "Dielectric sphere" begin
        # define scatterer: dielectric sphere
        sp = DielectricSphere(; radius=1.0, embedding=Medium(ε0, μ0), filling=Medium(ε0 * 5, μ0))

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
        sp = PECSphere(1.0, Medium(ε0, μ0))

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
        sp = LayeredSphere(; radii=SVector(1.0, 0.5), filling=SVector(Medium(3ε0, μ0), Medium(5ε0, μ0)), embedding=Medium(ε0, μ0))

        # define observation points in both layers and outside of the sphere
        point_cart = [SVector(0.25, 0.0, 0.0), SVector(0.75, 0.0, 0.0), SVector(2.0, 0.0, 0.0)]

        # compute scattered field and potential
        E = scatteredfield(sp, ex, ElectricField(point_cart))
        Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

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
        sp = LayeredSpherePEC(; radii=SVector(1.0, 0.5), embedding=Medium(ε0, μ0), filling=SVector(Medium(5ε0, μ0)))

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
    end
end
