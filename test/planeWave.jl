
@testset "PEC" begin

    f = 1e8

    sp = PECSphere(; radius=spRadius, embedding=Medium(𝜀, 𝜇))
    ex = planeWave(sp; frequency=f)

    @testset "Planewave excitation" begin
        @test planeWave(sp; frequency=f) isa PlaneWave{Float64,Float64,Float64}
    end

    @testset "Incident fields" begin

        point_cart = [SVector(2.0, 2.0, 3.2)]

        @test_nowarn E = field(ex, ElectricField(point_cart))
        @test_nowarn H = field(ex, MagneticField(point_cart))

    end

    @testset "Scattered fields" begin

        # ----- BEAST solution
        κ = 2π * f / c   # Wavenumber

        𝐸 = Maxwell3D.planewave(; direction=ẑ, polarization=x̂, wavenumber=κ)

        𝑒 = n × 𝐸 × n
        𝑇 = Maxwell3D.singlelayer(; wavenumber=κ, alpha=-im * 𝜇 * (2π * f), beta=1 / (-im * 𝜀 * (2π * f)))

        e = -assemble(𝑒, RT)
        T = assemble(𝑇, RT, RT)

        u = T \ e

        EF_MoM = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
        HF_MoM = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
        FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)

        # ----- this package
        EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
        HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
        FF = scatteredfield(sp, ex, FarField(points_cartFF))


        # ----- compare
        diff_EF = norm.(EF - EF_MoM) ./ maximum(norm.(EF))  # worst case error
        diff_HF = norm.(HF - HF_MoM) ./ maximum(norm.(HF))  # worst case error
        diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

        #@show maximum(20 * log10.(abs.(diff_EF)))
        #@show maximum(20 * log10.(abs.(diff_HF)))
        #@show maximum(20 * log10.(abs.(diff_FF)))

        @test maximum(20 * log10.(abs.(diff_EF))) < -25 # dB
        @test maximum(20 * log10.(abs.(diff_HF))) < -25 # dB
        @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
    end
end
