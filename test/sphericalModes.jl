
f = 1e8
κ = 2π * f / c   # Wavenumber

# BEAST impedance matrix
𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)
T = assemble(𝑇, RT, RT)


@testset "TE modes" begin

    @testset "Incident fields" begin

        ex = SphericalModeTE(; frequency=f, m=0, n=1, c=2)

        # define an observation point
        point_cart = [SVector(0.0, 1.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        # rather poor test: values are checked against values determined by this package
        #                   only to see whether future changes result in the same values
        @test FF[1][1] ≈ 0.005309361037222248 + 0.0im
        @test FF[1][2] ≈ -3.2510459998759485e-19 + 0.0im
        @test FF[1][3] ≈ 0.0 + 0.0im

        @test EF[1][1] ≈ 0.0010133416507273605 - 0.0012376919299125465im
        @test EF[1][2] ≈ -6.204928045029786e-20 + 7.578677301489551e-20im
        @test EF[1][3] ≈ 0.0 + 0.0im

        @test HF[1][1] ≈ 9.889059303803533e-23 - 2.333915413016954e-22im
        @test HF[1][2] ≈ 1.6150059446835905e-6 - 3.811573124009171e-6im
        @test HF[1][3] ≈ -3.6397996407963353e-6 - 1.3757077858915473e-6im


        ex3 = SphericalModeTE(; frequency=f, m=0, n=1, c=3)

        @test_throws ErrorException("Type can only be 1 or 2.") EF = field(ex3, ElectricField(point_cart))
        @test_throws ErrorException("Type can only be 1 or 2.") HF = field(ex3, MagneticField(point_cart))
    end

    @testset "Scattered fields" begin

        # check all different combinations of n and m till n = 4
        for nn in 1:4
            for mm in (-nn):nn

                𝐸 = ex = SphericalModeTE(; frequency=f, m=mm, n=nn, c=1)

                # ----- BEAST solution
                𝑒 = n × 𝐸 × n
                e = -assemble(𝑒, RT)

                u = T \ e

                EF_MoM = potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
                HF_MoM = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT) / (c * 𝜇)
                FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)

                # ----- this package
                sp = PECSphere(; radius=spRadius)

                EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
                HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
                FF = scatteredfield(sp, ex, FarField(points_cartFF))

                # ----- compare
                diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
                diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
                diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

                @test maximum(20 * log10.(abs.(diff_EF))) < -20 # dB 
                @test maximum(20 * log10.(abs.(diff_HF))) < -20 # dB
                @test maximum(20 * log10.(abs.(diff_FF))) < -20 # dB
            end
        end
    end

    @testset "Total fields" begin

        ex = SphericalModeTE(; frequency=f, m=0, n=1, c=2)
        sp = PECSphere(; radius=spRadius)

        # define an observation point
        point_cart = [SVector(2.0, 2.0, 3.2), SVector(3.1, 4, 2)]

        # compute scattered fields
        Es = scatteredfield(sp, ex, ElectricField(point_cart))
        Hs = scatteredfield(sp, ex, MagneticField(point_cart))

        Ei = field(ex, ElectricField(point_cart))
        Hi = field(ex, MagneticField(point_cart))

        # total field
        E = field(sp, ex, ElectricField(point_cart))
        H = field(sp, ex, MagneticField(point_cart))

        @test_throws ErrorException("The total far-field for a spherical mode excitation is not defined") field(
            sp, ex, FarField(point_cart)
        )

        # is it the sum?
        @test E[1] == Es[1] .+ Ei[1]
        @test H[1] == Hs[1] .+ Hi[1]
    end
end


@testset "TM modes" begin

    @testset "Incident fields" begin

        ex = SphericalModeTM(; frequency=f, m=0, n=1, c=2)

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        # rather poor test: values are checked against values determined by this package
        #                   only to see whether future changes result in the same values
        @test FF[1][1] ≈ 0.0 - 0.007534485420840261im
        @test FF[1][2] ≈ 0.0 - 0.0037672427104201303im
        @test FF[1][3] ≈ 0.0 + 0.011772633470062908im

        @test EF[1][1] ≈ 0.0009754932212024964 - 0.0009843946828623155im
        @test EF[1][2] ≈ 0.0004877466106012482 - 0.0004921973414311578im
        @test EF[1][3] ≈ -0.0018500711620225711 + 0.0010779970032754314im

        @test HF[1][1] ≈ -2.5597609870942054e-6 + 1.8128574477362876e-6im
        @test HF[1][2] ≈ 5.119521974188411e-6 - 3.625714895472575e-6im
        @test HF[1][3] ≈ 0.0 + 0.0im

        ex3 = SphericalModeTM(; frequency=f, m=0, n=1, c=3)

        @test_throws ErrorException("Type can only be 1 or 2.") EF = field(ex3, ElectricField(point_cart))
        @test_throws ErrorException("Type can only be 1 or 2.") HF = field(ex3, MagneticField(point_cart))
    end

    @testset "Scattered fields" begin

        # check all different combinations of n and m till n = 4
        for nn in 1:4
            for mm in (-nn):nn

                𝐸 = ex = SphericalModeTM(; frequency=f, m=mm, n=nn, c=1)

                # ----- BEAST solution
                𝑒 = n × 𝐸 × n
                e = -assemble(𝑒, RT)

                u = T \ e

                EF_MoM = potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
                HF_MoM = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT) / (c * 𝜇)
                FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(; gamma=𝑇.gamma), points_cartFF, u, RT)

                # ----- this package
                sp = PECSphere(; radius=spRadius)

                EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
                HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
                FF = scatteredfield(sp, ex, FarField(points_cartFF))

                # ----- compare
                diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
                diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
                diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

                @test maximum(20 * log10.(abs.(diff_EF))) < -20 # dB 
                @test maximum(20 * log10.(abs.(diff_HF))) < -20 # dB
                @test maximum(20 * log10.(abs.(diff_FF))) < -20 # dB
            end
        end
    end

    @testset "Total fields" begin

        ex = SphericalModeTM(; frequency=f, m=0, n=1, c=2)
        sp = PECSphere(; radius=spRadius)

        # define an observation point
        point_cart = [SVector(2.0, 2.0, 3.2), SVector(3.1, 4, 2)]

        # compute scattered fields
        Es = scatteredfield(sp, ex, ElectricField(point_cart))
        Hs = scatteredfield(sp, ex, MagneticField(point_cart))

        Ei = field(ex, ElectricField(point_cart))
        Hi = field(ex, MagneticField(point_cart))

        # total field
        E = field(sp, ex, ElectricField(point_cart))
        H = field(sp, ex, MagneticField(point_cart))

        @test_throws ErrorException("The total far-field for a spherical mode excitation is not defined") field(
            sp, ex, FarField(point_cart)
        )

        # is it the sum?
        @test E[1] == Es[1] .+ Ei[1]
        @test H[1] == Hs[1] .+ Hi[1]
    end
end
