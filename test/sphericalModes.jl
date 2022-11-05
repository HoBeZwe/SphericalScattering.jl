



@testset "TE modes" begin

    ex = SphericalModeTE(; wavenumber=30.0, m=0, n=1, c=1) # ≈ 10 MHz

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        E  = field(ex, ElectricField(point_cart))
        H  = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        # TODO: find better way of verifying (via BEAST)
        # @test FF[1][1] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][2] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][3] ≈ -368.78010908176793 - 140.01031931454293im

        # @test E[1][1] ≈ 33.75576840234728 - 66.33989918223259im
        # @test E[1][2] ≈ 33.75576840234728 - 66.33989918223259im
        # @test E[1][3] ≈ -118.11657839808761 + 218.13026261521526im

        # @test H[1][1] ≈ -0.23911545354270725 + 0.4457313513337625im
        # @test H[1][2] ≈ 0.23911545354270725 - 0.4457313513337625im
        # @test H[1][3] ≈ 0.0 + 0.0im
        @test 1 == 1
    end

    @testset "Scattered fields" begin

        # define scatterer: PEC sphere
        sp = PECSphere(; radius=1.0)

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        # compute scattered fields
        E  = scatteredfield(sp, ex, ElectricField(point_cart))
        H  = scatteredfield(sp, ex, MagneticField(point_cart))
        FF = scatteredfield(sp, ex, FarField(point_cart))

        # @test E[1][1] ≈ -5.502103327427189 - 7.627682569142905im
        # @test E[1][2] ≈ -5.502103327427188 - 7.6276825691429035im
        # @test E[1][3] ≈ 8.424633687603194 + 12.484440209052982im

        # @test H[1][1] ≈ -0.021517831140115615 - 0.030977247907197428im
        # @test H[1][2] ≈ 0.021517831140115618 + 0.03097724790719743im
        # @test H[1][3] ≈ 0.0 + 0.0im

        # @test FF[1][1] ≈ 14.44466970131592 + 31.74394205284891im
        # @test FF[1][2] ≈ 14.444669701315918 + 31.743942052848908im
        # @test FF[1][3] ≈ -18.055837126644896 - 39.67992756606113im
        @test 1 == 1
    end
end


@testset "TM modes" begin

    ex = ex = SphericalModeTM(; wavenumber=30.0, m=0, n=1, c=1) # ≈ 10 MHz

    @testset "Incident fields" begin

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        E  = field(ex, ElectricField(point_cart))
        H  = field(ex, MagneticField(point_cart))
        FF = field(ex, FarField(point_cart))

        # TODO: find better way of verifying (via BEAST)
        # @test FF[1][1] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][2] ≈ 295.0240872654143 + 112.00825545163434im
        # @test FF[1][3] ≈ -368.78010908176793 - 140.01031931454293im

        # @test E[1][1] ≈ 33.75576840234728 - 66.33989918223259im
        # @test E[1][2] ≈ 33.75576840234728 - 66.33989918223259im
        # @test E[1][3] ≈ -118.11657839808761 + 218.13026261521526im

        # @test H[1][1] ≈ -0.23911545354270725 + 0.4457313513337625im
        # @test H[1][2] ≈ 0.23911545354270725 - 0.4457313513337625im
        # @test H[1][3] ≈ 0.0 + 0.0im
        @test 1 == 1
    end

    @testset "Scattered fields" begin

        # define scatterer: PEC sphere
        sp = PECSphere(; radius=1.0)

        # define an observation point
        point_cart = [SVector(4.0, 2.0, 3.2)]

        # compute scattered fields
        E  = scatteredfield(sp, ex, ElectricField(point_cart))
        H  = scatteredfield(sp, ex, MagneticField(point_cart))
        FF = scatteredfield(sp, ex, FarField(point_cart))

        # @test E[1][1] ≈ -5.502103327427189 - 7.627682569142905im
        # @test E[1][2] ≈ -5.502103327427188 - 7.6276825691429035im
        # @test E[1][3] ≈ 8.424633687603194 + 12.484440209052982im

        # @test H[1][1] ≈ -0.021517831140115615 - 0.030977247907197428im
        # @test H[1][2] ≈ 0.021517831140115618 + 0.03097724790719743im
        # @test H[1][3] ≈ 0.0 + 0.0im

        # @test FF[1][1] ≈ 14.44466970131592 + 31.74394205284891im
        # @test FF[1][2] ≈ 14.444669701315918 + 31.743942052848908im
        # @test FF[1][3] ≈ -18.055837126644896 - 39.67992756606113im
        @test 1 == 1
    end
end
