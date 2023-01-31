
@testset "PEC" begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    sp = PECSphere(; radius=spRadius, embedding=Medium(𝜀, 𝜇))
    ex = planeWave(sp; frequency=f)

    @testset "Planewave excitation" begin
        @test planeWave(; frequency=f) isa PlaneWave{Float64,Float64,Float64}
        @test planeWave(sp; frequency=f) isa PlaneWave{Float64,Float64,Float64}
    end

    @testset "Incident fields" begin

        point_cart = [SVector(2.0, 2.0, 3.2)]

        @test_nowarn E = field(ex, ElectricField(point_cart))
        @test_nowarn H = field(ex, MagneticField(point_cart))

        @test_throws ErrorException("The far-field of a plane wave is not defined.") field(ex, FarField(point_cart))

    end

    @testset "Scattered fields" begin

        @testset "Standard orientation" begin

            # ----- BEAST solution
            𝐸 = Maxwell3D.planewave(; direction=ẑ, polarization=x̂, wavenumber=κ)

            𝑒 = n × 𝐸 × n
            𝑇 = Maxwell3D.singlelayer(; wavenumber=κ, alpha=-im * 𝜇 * (2π * f), beta=1 / (-im * 𝜀 * (2π * f)))

            e = -assemble(𝑒, RT)
            T = assemble(𝑇, RT, RT)

            u = T \ e

            EF_MoM₂ = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
            HF_MoM₂ = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
            FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)

            # ----- this package
            ex = planeWave(sp; frequency=f)

            EF₂ = scatteredfield(sp, ex, ElectricField(points_cartNF))
            EF₁ = scatteredfield(sp, ex, ElectricField(points_cartNF_inside))
            HF₂ = scatteredfield(sp, ex, MagneticField(points_cartNF))
            HF₁ = scatteredfield(sp, ex, MagneticField(points_cartNF_inside))
            FF = scatteredfield(sp, ex, FarField(points_cartFF))


            # ----- compare
            diff_EF₂ = norm.(EF₂ - EF_MoM₂) ./ maximum(norm.(EF₂))  # worst case error
            diff_HF₂ = norm.(HF₂ - HF_MoM₂) ./ maximum(norm.(HF₂))  # worst case error
            diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

            @test maximum(20 * log10.(abs.(diff_EF₂))) < -25 # dB
            @test norm(EF₁) == 0.0
            @test maximum(20 * log10.(abs.(diff_HF₂))) < -25 # dB
            @test norm(HF₁) == 0.0
            @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
        end

        @testset "General orientation" begin

            # ----- BEAST solution
            dir = normalize(SVector(0.0, 1.0, 1.0)) # normalization for BEAST
            pol = normalize(SVector(-1.0, 0.0, 0.0))


            𝐸 = Maxwell3D.planewave(; direction=dir, polarization=pol, wavenumber=κ)

            𝑒 = n × 𝐸 × n
            𝑇 = Maxwell3D.singlelayer(; wavenumber=κ, alpha=-im * 𝜇 * (2π * f), beta=1 / (-im * 𝜀 * (2π * f)))

            e = -assemble(𝑒, RT)
            T = assemble(𝑇, RT, RT)

            u = T \ e

            EF_MoM₂ = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
            HF_MoM₂ = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
            FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)

            # ----- this package
            ex = planeWave(sp; frequency=f, direction=dir, polarization=pol)

            EF₂ = scatteredfield(sp, ex, ElectricField(points_cartNF))
            EF₁ = scatteredfield(sp, ex, ElectricField(points_cartNF_inside))
            HF₂ = scatteredfield(sp, ex, MagneticField(points_cartNF))
            HF₁ = scatteredfield(sp, ex, MagneticField(points_cartNF_inside))
            FF = scatteredfield(sp, ex, FarField(points_cartFF))


            # ----- compare
            diff_EF₂ = norm.(EF₂ - EF_MoM₂) ./ maximum(norm.(EF₂))  # worst case error
            diff_HF₂ = norm.(HF₂ - HF_MoM₂) ./ maximum(norm.(HF₂))  # worst case error
            diff_FF = norm.(FF - FF_MoM) ./ maximum(norm.(FF))  # worst case error

            @test maximum(20 * log10.(abs.(diff_EF₂))) < -25 # dB
            @test norm(EF₁) == 0.0
            @test maximum(20 * log10.(abs.(diff_HF₂))) < -25 # dB
            @test norm(HF₁) == 0.0
            @test maximum(20 * log10.(abs.(diff_FF))) < -25 # dB
        end
    end


    @testset "Total fields" begin

        # define an observation point
        point_cart = [SVector(2.0, 2.0, 3.2), SVector(3.1, 4, 2)]

        # compute scattered fields
        Es = scatteredfield(sp, ex, ElectricField(point_cart))
        Hs = scatteredfield(sp, ex, MagneticField(point_cart))
        #FFs = scatteredfield(sp, ex, FarField(point_cart))

        Ei = field(ex, ElectricField(point_cart))
        Hi = field(ex, MagneticField(point_cart))
        #FFi = field(ex, FarField(point_cart))

        # total field
        E = field(sp, ex, ElectricField(point_cart))
        H = field(sp, ex, MagneticField(point_cart))
        @test_throws ErrorException("The total far-field for a plane-wave excitation is not defined") field(
            sp, ex, FarField(point_cart)
        )

        # is it the sum?
        @test E[1] == Es[1] .+ Ei[1]
        @test H[1] == Hs[1] .+ Hi[1]
    end
end
