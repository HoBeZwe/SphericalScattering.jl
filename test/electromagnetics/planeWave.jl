
@testitem "PEC" setup = [Setup, BEASTSetup] begin

    f = 1e8
    κ = 2π * f / c   # Wavenumber

    sp = PECSphere(; radius=spRadius)
    ex = planeWave(; frequency=f)

    @testset "Planewave excitation" begin
        @test planeWave(; frequency=f) isa PlaneWave{Float64,Float64,Float64}
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
            𝑇 = Maxwell3D.singlelayer(; wavenumber=κ, alpha=(-im * 𝜇 * (2π * f)), beta=1 / (-im * 𝜀 * (2π * f)))

            e = -assemble(𝑒, RT)
            T = assemble(𝑇, RT, RT)

            u = T \ e

            EF_MoM₂ = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
            HF_MoM₂ = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
            FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)

            # ----- this package
            ex = planeWave(; frequency=f)

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
            𝑇 = Maxwell3D.singlelayer(; wavenumber=κ, alpha=(-im * 𝜇 * (2π * f)), beta=1 / (-im * 𝜀 * (2π * f)))

            e = -assemble(𝑒, RT)
            T = assemble(𝑇, RT, RT)

            u = T \ e

            EF_MoM₂ = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
            HF_MoM₂ = potential(BEAST.MWDoubleLayerField3D(; wavenumber=κ), points_cartNF, u, RT)
            FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)

            # ----- this package
            ex = planeWave(; frequency=f, direction=dir, polarization=pol)

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


@testitem "PEC boundary conditions and limits" setup = [Setup] begin

    @testset "Boundary conditions" begin

        # On the surface of a PEC sphere the total field has to satisfy n × E = 0 and n ⋅ H = 0.
        # The limit from the outside is approximated at the (relative) distance δ from the surface,
        # which leaves an error in the order of (1 + ka) * δ.
        δ = 1e-12
        tol = 1e-9

        function boundaryErrors(sp, ex)

            points_cart, ~ = getDefaultPoints(sp.radius * (1 + δ))
            normals = normalize.(points_cart)

            # ----- incident fields
            Ei = field(ex, ElectricField(points_cart))
            Hi = field(ex, MagneticField(points_cart))

            # ----- total fields
            E = field(sp, ex, ElectricField(points_cart))
            H = field(sp, ex, MagneticField(points_cart))

            diff_Et = norm.(cross.(normals, E)) ./ maximum(norm.(Ei))
            diff_Hn = abs.(dot.(normals, H)) ./ maximum(norm.(Hi))

            return maximum(diff_Et), maximum(diff_Hn)
        end

        function insideFields(sp, ex)

            points_cart, ~ = getDefaultPoints(sp.radius / 2)

            E = field(sp, ex, ElectricField(points_cart))
            H = field(sp, ex, MagneticField(points_cart))

            return E, H
        end

        # electrical sizes from ka ≈ 0.002 to ka ≈ 63
        @testset "Electrical size: f = $f Hz, radius = $radius m" for f in (1e6, 1e8, 1e9), radius in (0.1, 1.0, 3.0)

            sp = PECSphere(; radius=radius)
            ex = planeWave(; frequency=f)

            diff_Et, diff_Hn = boundaryErrors(sp, ex)
            E₁, H₁ = insideFields(sp, ex)

            @test diff_Et < tol
            @test diff_Hn < tol
            @test norm(E₁) == 0.0
            @test norm(H₁) == 0.0
        end

        @testset "General orientation" begin

            orientations = (
                (SVector(0.0, 1.0, 1.0), SVector(-1.0, 0.0, 0.0)),
                (SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)),
                (SVector(0.0, 0.0, -1.0), SVector(0.0, 1.0, 0.0)),
            )

            for (dir, pol) in orientations

                sp = PECSphere(; radius=spRadius)
                ex = planeWave(; frequency=f, direction=dir, polarization=pol)

                diff_Et, diff_Hn = boundaryErrors(sp, ex)
                E₁, H₁ = insideFields(sp, ex)

                @test diff_Et < tol
                @test diff_Hn < tol
                @test norm(E₁) == 0.0
                @test norm(H₁) == 0.0
            end
        end

        @testset "Embedding and amplitude" begin

            sp = PECSphere(; radius=spRadius)
            ex = planeWave(; frequency=f, embedding=Medium(𝜀 * 3.0, 𝜇 * 2.0), amplitude=2.5)

            diff_Et, diff_Hn = boundaryErrors(sp, ex)
            E₁, H₁ = insideFields(sp, ex)

            @test diff_Et < tol
            @test diff_Hn < tol
            @test norm(E₁) == 0.0
            @test norm(H₁) == 0.0
        end
    end

    relativeDifference(F, Fref) = maximum(norm.(F - Fref)) / maximum(norm.(Fref))

    @testset "RCS limits" begin

        radius = 2.0
        sp = PECSphere(; radius=radius)

        frequencyOf(ka) = ka * c / (2π * radius)

        # Electrically small sphere, the monostatic RCS is 9πa² (ka)⁴
        @testset "Rayleigh limit: ka = $ka" for ka in (0.001, 0.01, 0.05)
            σ = rcs(sp, planeWave(; frequency=frequencyOf(ka)))

            @test isapprox(σ, 9π * radius^2 * ka^4; rtol=ka^2)
        end

        # Electrically large sphere, the monostatic RCS approaches the geometrical cross section πa².
        @testset "Optical limit: ka = $ka" for ka in (100.0, 200.0, 500.0)
            σ = rcs(sp, planeWave(; frequency=frequencyOf(ka)))

            @test isapprox(σ, π * radius^2; rtol=1e-2)
        end
    end

    @testset "Rotation of the plane wave" begin

        sp = PECSphere(; radius=spRadius)
        ex₀ = planeWave(; frequency=f)

        E₀ = field(sp, ex₀, ElectricField(points_cartNF))
        H₀ = field(sp, ex₀, MagneticField(points_cartNF))
        FF₀ = scatteredfield(sp, ex₀, FarField(points_cartFF))

        orientations = (
            (SVector(0.0, 1.0, 1.0), SVector(-1.0, 0.0, 0.0)),
            (SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)),
            (SVector(0.0, 0.0, -1.0), SVector(0.0, 1.0, 0.0)),
            (SVector(1.0, 2.0, 2.0), SVector(2.0, 1.0, -2.0)),
        )

        for (dir, pol) in orientations

            ex = planeWave(; frequency=f, direction=dir, polarization=pol)

            d = normalize(dir)
            p = normalize(pol)
            R = hcat(p, cross(d, p), d) # maps x̂ → p, ŷ → d × p, ẑ → d

            rotated(vectors) = [R * vector for vector in vectors]

            E = field(sp, ex, ElectricField(rotated(points_cartNF)))
            H = field(sp, ex, MagneticField(rotated(points_cartNF)))
            FF = scatteredfield(sp, ex, FarField(rotated(points_cartFF)))

            @test relativeDifference(E, rotated(E₀)) < 1e-12
            @test relativeDifference(H, rotated(H₀)) < 1e-12
            @test relativeDifference(FF, rotated(FF₀)) < 1e-12
            @test rcs(sp, ex) ≈ rcs(sp, ex₀)
        end
    end

    @testset "Linearity in the amplitude" begin

        # The fields scale with the amplitude of the plane wave, the RCS does not depend on it.
        sp = PECSphere(; radius=spRadius)
        ex₁ = planeWave(; frequency=f)

        E₁ = field(sp, ex₁, ElectricField(points_cartNF))
        H₁ = field(sp, ex₁, MagneticField(points_cartNF))
        FF₁ = scatteredfield(sp, ex₁, FarField(points_cartFF))

        for amplitude in (1e-3, 2.5, 40.0)

            ex = planeWave(; frequency=f, amplitude=amplitude)

            E = field(sp, ex, ElectricField(points_cartNF))
            H = field(sp, ex, MagneticField(points_cartNF))
            FF = scatteredfield(sp, ex, FarField(points_cartFF))

            @test relativeDifference(E, amplitude * E₁) < 1e-12
            @test relativeDifference(H, amplitude * H₁) < 1e-12
            @test relativeDifference(FF, amplitude * FF₁) < 1e-12
            @test rcs(sp, ex) ≈ rcs(sp, ex₁)
        end
    end
end
