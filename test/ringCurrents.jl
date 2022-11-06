

@testset "Electric ring current" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber


    ex = electricRingCurrent(; wavenumber=κ, center=SVector(0.0, 0.0, 2.0), radius=0.5)

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        @test FF[1][1] ≈ 7.942044920729696 - 3.1030988210509873im
        @test FF[1][2] ≈ -15.88408984145939 + 6.206197642101974im
        @test FF[1][3] ≈ 0.0 + 0.0im

        @test EF[1][1] ≈ -7.942044920729696 + 3.1030988210509873im
        @test EF[1][2] ≈ 15.88408984145939 - 6.206197642101974im
        @test EF[1][3] ≈ 0.0 + 0.0im

        @test HF[1][1] ≈ -0.01019515697859331 + 0.00628333403389956im
        @test HF[1][2] ≈ -0.0050975784892966555 + 0.003141667016949781im
        @test HF[1][3] ≈ 0.04525565726450052 - 0.017077464855190773im
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
        HF_MoM = -potential(MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)


        # ----- this package
        sp = PECSphere(; radius=spRadius)

        EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
        HF = scatteredfield(sp, ex, MagneticField(points_cartNF)) * c * 𝜇
        FF = scatteredfield(sp, ex, FarField(points_cartFF))


        # ----- compare
        diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
        diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
        diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

        @test maximum(20 * log10.(abs.(diff_EF))) < -24 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -24 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -24 # dB
    end
end


@testset "Magnetic ring current" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber


    ex = magneticRingCurrent(; wavenumber=κ, center=SVector(0.0, 0.0, 2.0), radius=0.5)

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        @test FF[1][1] ≈ -5.595916238822697e-5 + 2.186424435609382e-5im
        @test FF[1][2] ≈ 0.00011191832477645392 - 4.3728488712187633e-5im
        @test FF[1][3] ≈ 0.0 + 0.0im

        @test EF[1][1] ≈ 0.01019515697859331 - 0.00628333403389956im
        @test EF[1][2] ≈ 0.0050975784892966555 - 0.003141667016949781im
        @test EF[1][3] ≈ -0.04525565726450052 + 0.017077464855190773im

        @test HF[1][1] ≈ -5.595916238822697e-5 + 2.186424435609382e-5im
        @test HF[1][2] ≈ 0.00011191832477645392 - 4.3728488712187633e-5im
        @test HF[1][3] ≈ 0.0 + 0.0im
    end

    @testset "Scattered fields" begin

        # ----- BEAST solution
        𝐸 = ex

        𝑒 = n × 𝐸 × n
        𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)

        e = assemble(𝑒, RT)
        T = assemble(𝑇, RT, RT)

        u = T \ e

        EF_MoM = -potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        HF_MoM = +potential(MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        FF_MoM = +im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)


        # ----- this package
        sp = PECSphere(; radius=spRadius)

        EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
        HF = scatteredfield(sp, ex, MagneticField(points_cartNF)) * c * 𝜇
        FF = scatteredfield(sp, ex, FarField(points_cartFF))


        # ----- compare
        diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
        diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
        diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

        @test maximum(20 * log10.(abs.(diff_EF))) < -24 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -24 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -24 # dB
    end
end
