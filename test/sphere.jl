
@testset "Medium and Sphere" begin

    md = Medium(3.0, Float32(2.0))

    @test md isa Medium{Float64}

    mdf = Medium(Float32(2), Float32(1))

    @test mdf isa Medium{Float32}

    @test Medium{Float64}(mdf) isa Medium{Float64}

    sp = DielectricSphere(; radius=Float32(1.0), embedding=mdf, filling=md)

    @test sp isa DielectricSphere{Float64,Float32}

end
