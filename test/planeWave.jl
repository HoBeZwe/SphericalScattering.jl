
@testset "Plane wave" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    ex = planeWave(; wavenumber=κ)


    @testset "Incident fields" begin

        point_cart = [SVector(2.0, 2.0, 3.2)]

        @test_nowarn E = field(ex, ElectricField(point_cart))
        @test_nowarn H = field(ex, MagneticField(point_cart))

    end

    @testset "Scattered fields" begin

        # ----- BEAST solution
        𝐸 = Maxwell3D.planewave(; direction=-ẑ, polarization=x̂, wavenumber=κ)

        𝑒 = n × 𝐸 × n
        𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)

        e = assemble(𝑒, RT)
        T = assemble(𝑇, RT, RT)

        u = T \ e

        EF_MoM = potential(MWSingleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        HF_MoM = potential(MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
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

        #@show maximum(20 * log10.(abs.(diff_EF)))
        #@show maximum(20 * log10.(abs.(diff_HF)))
        #@show maximum(20 * log10.(abs.(diff_FF)))

        @test maximum(20 * log10.(abs.(diff_EF))) < -27 # dB 
        @test maximum(20 * log10.(abs.(diff_HF))) < -27 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -27 # dB
    end
end
