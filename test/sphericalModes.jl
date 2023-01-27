
@testset "TE modes" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber


    ex = SphericalModeTE(; frequency=f, m=0, n=1, c=2)

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        @test_nowarn EF = field(ex, ElectricField(point_cart))
        @test_nowarn HF = field(ex, MagneticField(point_cart))
        @test_nowarn FF = field(ex, FarField(point_cart))

        point_cart = [SVector(0.0, 0.0, 3.2), SVector(0.0, 0.0, -3.2)]

        ex2 = SphericalModeTE(; frequency=f, m=1, n=1, c=1)

        @test_nowarn EF = field(ex2, ElectricField(point_cart))
        @test_nowarn HF = field(ex2, MagneticField(point_cart))
        @test_nowarn FF = field(ex2, FarField(point_cart))

        ex3 = SphericalModeTE(; frequency=f, m=0, n=1, c=3)

        @test_throws ErrorException("Type can only be 1 or 2.") EF = field(ex3, ElectricField(point_cart))
        @test_throws ErrorException("Type can only be 1 or 2.") HF = field(ex3, MagneticField(point_cart))

        # @test FF[1][1] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][2] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][3] ≈ -368.78010908176793 - 140.01031931454293im

        # @test EF[1][1] ≈ 33.75576840234728 - 66.33989918223259im
        # @test EF[1][2] ≈ 33.75576840234728 - 66.33989918223259im
        # @test EF[1][3] ≈ -118.11657839808761 + 218.13026261521526im

        # @test H[1][1] ≈ -0.23911545354270725 + 0.4457313513337625im
        # @test H[1][2] ≈ 0.23911545354270725 - 0.4457313513337625im
        # @test H[1][3] ≈ 0.0 + 0.0im
    end

    @testset "Scattered fields" begin

        # ----- BEAST solution
        𝐸 = ex

        𝑒 = n × 𝐸 × n
        𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)

        e = assemble(𝑒, RT)
        T = assemble(𝑇, RT, RT)

        u = T \ e

        EF_MoM = +potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        HF_MoM = -potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)


        # ----- this package
        sp = PECSphere(; radius=spRadius)

        EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
        HF = scatteredfield(sp, ex, MagneticField(points_cartNF)) * c * 𝜇
        FF = scatteredfield(sp, ex, FarField(points_cartFF))

        # ----- compare
        #diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
        #diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
        #diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

        diff_EF = norm.(EF) ./ maximum(norm.(EF)) - norm.(EF_MoM) ./ maximum(norm.(EF_MoM))  # worst case error
        diff_HF = norm.(HF) ./ maximum(norm.(HF)) - norm.(HF_MoM) ./ maximum(norm.(HF_MoM))  # worst case error
        diff_FF = norm.(FF) ./ maximum(norm.(FF)) - norm.(FF_MoM) ./ maximum(norm.(FF_MoM))  # worst case error

        @test maximum(20 * log10.(abs.(diff_EF))) < -25 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -25 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
    end
end


@testset "TM modes" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber


    ex = SphericalModeTM(; frequency=f, m=0, n=1, c=2)

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        @test_nowarn EF = field(ex, ElectricField(point_cart))
        @test_nowarn HF = field(ex, MagneticField(point_cart))
        @test_nowarn FF = field(ex, FarField(point_cart))

        ex2 = SphericalModeTM(; frequency=f, m=1, n=1, c=1)

        point_cart = [SVector(0.0, 0.0, 3.2), SVector(0.0, 0.0, -3.2)]

        @test_nowarn EF = field(ex2, ElectricField(point_cart))
        @test_nowarn HF = field(ex2, MagneticField(point_cart))
        @test_nowarn FF = field(ex2, FarField(point_cart))

        ex3 = SphericalModeTM(; frequency=f, m=0, n=1, c=3)

        @test_throws ErrorException("Type can only be 1 or 2.") EF = field(ex3, ElectricField(point_cart))
        @test_throws ErrorException("Type can only be 1 or 2.") HF = field(ex3, MagneticField(point_cart))

        # @test FF[1][1] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][2] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][3] ≈ -368.78010908176793 - 140.01031931454293im

        # @test EF[1][1] ≈ 33.75576840234728 - 66.33989918223259im
        # @test EF[1][2] ≈ 33.75576840234728 - 66.33989918223259im
        # @test EF[1][3] ≈ -118.11657839808761 + 218.13026261521526im

        # @test HF[1][1] ≈ -0.23911545354270725 + 0.4457313513337625im
        # @test HF[1][2] ≈ 0.23911545354270725 - 0.4457313513337625im
        # @test HF[1][3] ≈ 0.0 + 0.0im
    end

    @testset "Scattered fields" begin

        # ----- BEAST solution
        𝐸 = ex

        𝑒 = n × 𝐸 × n
        𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)

        e = assemble(𝑒, RT)
        T = assemble(𝑇, RT, RT)

        u = T \ e

        EF_MoM = +potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        HF_MoM = -potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)


        # ----- this package
        sp = PECSphere(; radius=spRadius)

        EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
        HF = scatteredfield(sp, ex, MagneticField(points_cartNF)) * c * 𝜇
        FF = scatteredfield(sp, ex, FarField(points_cartFF))


        # ----- compare
        #diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
        #diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
        #diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

        diff_EF = norm.(EF) ./ maximum(norm.(EF)) - norm.(EF_MoM) ./ maximum(norm.(EF_MoM))  # worst case error
        diff_HF = norm.(HF) ./ maximum(norm.(HF)) - norm.(HF_MoM) ./ maximum(norm.(HF_MoM))  # worst case error
        diff_FF = norm.(FF) ./ maximum(norm.(FF)) - norm.(FF_MoM) ./ maximum(norm.(FF_MoM))  # worst case error

        @test maximum(20 * log10.(abs.(diff_EF))) < -25 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -25 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
    end
end
