
@testset "Hertzian dipole" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    ex = HertzianDipole(; frequency=f, center=SVector(0.0, 0.0, 2.0))


    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(2.0, 2.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        @test FF[1][1] ≈ -14.599166142315502 - 16.519749232363072im
        @test FF[1][2] ≈ -14.599166142315502 - 16.519749232363072im
        @test FF[1][3] ≈ 18.248957677894378 + 20.64968654045384im

        @test EF[1][1] ≈ -4.3874439814528525 + 3.1430494710254244im
        @test EF[1][2] ≈ -4.3874439814528525 + 3.1430494710254244im
        @test EF[1][3] ≈ 16.58730908095325 - 4.356087647969305im

        @test HF[1][1] ≈ 0.03405123034007954 - 0.010917855361887753im
        @test HF[1][2] ≈ -0.03405123034007954 + 0.010917855361887753im
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

        EF_MoM = +potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        HF_MoM = -potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
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

        @test maximum(20 * log10.(abs.(diff_EF))) < -25 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -25 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
    end
end


@testset "Fitzgerald dipole" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber


    ex = FitzgeraldDipole(; frequency=f, center=SVector(0.0, 0.0, 2.0))

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(2.0, 2.0, 3.2)]

        EF = field(ex, ElectricField(point_cart))
        HF = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        @test EF[1][1] ≈ -0.03405123034007954 + 0.010917855361887753im
        @test EF[1][2] ≈ 0.03405123034007954 - 0.010917855361887753im
        @test EF[1][3] ≈ 0.0 - 0.0im

        @test HF[1][1] ≈ -3.091366174302772e-5 + 2.2145734190481372e-5im
        @test HF[1][2] ≈ -3.091366174302772e-5 + 2.2145734190481372e-5im
        @test HF[1][3] ≈ 0.00011687316449470456 - 3.069272693658695e-5im

        @test FF[1][1] ≈ -0.0517201674997236 - 0.05852417795799779im
        @test FF[1][2] ≈ 0.0517201674997236 + 0.05852417795799779im
        @test FF[1][3] ≈ 0.0 + 0.0im
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
        HF_MoM = +potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
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
