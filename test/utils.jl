
@testset "Plotting points" begin

    # just check whether function calls returns no errors
    @test_nowarn sphericalGridPoints()
    @test_nowarn phiCutPoints(2)
    @test_nowarn thetaCutPoints(20)
end

@testset "Plotting extensions" begin

    using PlotlyJS

    # --- excite PEC sphere by Hertzian dipole
    orient = normalize(SVector(0.0, 1.0, 1.0))
    ex = HertzianDipole(; frequency=1e8, orientation=orient, position=2 * orient)
    sp = PECSphere(; radius=1.0)

    # --- evaluate fields at spherical grid points
    points_cart, points_sph = sphericalGridPoints()
    FF = scatteredfield(sp, ex, FarField(points_cart))

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=true, type="abs")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=true, type="abs")

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=true, type="theta")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=true, type="phi")

    @test_nowarn plotff(FF, points_sph, scale="linear", normalize=false, type="theta")
    @test_nowarn plotff(FF, points_sph, scale="log", normalize=false, type="phi")

    # --- evaluate fields at φ = 5° cut
    points_cart, points_sph = phiCutPoints(5) # analogously, thetaCutPoints can be used
    FF = scatteredfield(sp, ex, FarField(points_cart))

    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="polar")
    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="linear", format="polar")

    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="rectangular")
    @test_nowarn plotffcut(norm.(FF), points_sph, normalize=true, scale="linear", format="rectangular")
end

@testset "Radar cross section" begin


end
