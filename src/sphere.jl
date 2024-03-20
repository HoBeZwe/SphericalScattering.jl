"""
    Medium(ε, μ)

Homogeneous, isotropic background medium.
"""
struct Medium{C}
    ε::C
    μ::C
end

function Medium(ε::T1, μ::T2) where {T1,T2}
    T = promote_type(T1, T2)
    return Medium(T(ε), T(μ))
end

function Medium{T}(md) where {T}
    return Medium(T(md.ε), T(md.μ))
end

function Medium{T}(md::SVector{N,Medium{C}}) where {T,N,C}
    md2 = []
    for n in axes(md,1)
        println(md[n])
        push!(md2,Medium(T(md[n].ε),T(md[n].μ)))
    end
    return SVector{N,Medium{T}}(md2)
end 



abstract type Sphere end

struct DielectricSphere{C,R} <: Sphere
    radius::R
    filling::Medium{C}
end


"""
    DielectricSphere(
        radius      = error("missing argument `radius`"), 
        filling     = error("missing argument `filling`")
    )

Constructor for the dielectric sphere.
"""
DielectricSphere(; radius=error("missing argument `radius`"), filling=error("missing argument `filling`")) =
    DielectricSphere(radius, filling)



struct DielectricSphereThinImpedanceLayer{R,C} <: Sphere
    radius::R
    thickness::R
    thinlayer::Medium{C}
    filling::Medium{C}
end

function DielectricSphereThinImpedanceLayer(r::R1, d::R2, thinlayer::Medium{C1}, filling::Medium{C2}) where {R1,R2,C1,C2}

    R = promote_type(R1, R2)
    C = promote_type(C1, C2)

    DielectricSphereThinImpedanceLayer(R(r), R(d), Medium{C}(thinlayer), Medium{C}(filling))
end

"""
    DielectricSphereThinImpedanceLayer(
        radius    = error("missing argument `radius`"),
        thickness = error("missing argument `thickness` of the coating"),
        thinlayer = error("missing argument `thinlayer`"),
        filling   = error("missing argument `filling`")
    )

Constructor for the dielectric sphere with a thin impedance layer.
For this model, it is assumed that the displacement field is only radial
direction in the layer, which requires a small thickness and low conductivity.
For details, see for example T. B. Jones, Ed., “Models for layered spherical particles,”
in Electromechanics of Particles, Cambridge: Cambridge University Press, 1995, 
pp. 227–235. doi: 10.1017/CBO9780511574498.012.
"""
DielectricSphereThinImpedanceLayer(;
    radius=error("missing argument `radius`"),
    thickness=error("missing argument `thickness` of the coating"),
    thinlayer=error("missing argument `thinlayer`"),
    filling=error("missing argument `filling`"),
) = DielectricSphereThinImpedanceLayer(radius, thickness, thinlayer, filling)



struct PECSphere{R} <: Sphere
    radius::R
end

"""
    PECSphere( 
        radius = error("missing argument `radius`"), 
    )

Constructor for the PEC sphere.
"""
PECSphere(; radius=error("missing argument `radius`")) = PECSphere(radius)



struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{N,Medium{C}}
end

"""
    LayeredSphere( 
        radii   = error("Missing argument `radii`"), 
        filling = error("`missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
LayeredSphere(; radii=error("Missing argument `radii`"), filling=error("`missing argument `filling`")) = LayeredSphere(radii, filling)

function LayeredSphere(radii::SVector{N,R}, embedding::Medium{C1}, filling::SVector{N,Medium{C2}}) where {N,R,C1,C2}

    C = promote_type(C1, C2)

    LayeredSphere(radii, Medium{C}(embedding), Medium{C}(filling))
end


struct LayeredSpherePEC{N,D,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{D,Medium{C}}
end

"""
    LayeredSpherePEC( 
        radii   = error("Missing argument `radii`"), 
        filling = error("Missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
LayeredSpherePEC(; radii=error("Missing argument `radii`"), filling=error("Missing argument `filling`")) =
    LayeredSpherePEC(radii, filling)




"""
    layer(sp::Sphere, r)

Returns the index of the layer, `r` is located, where `1` denotes the inner most layer. 
"""
function layer(sp::Sphere, r)
    # Using Jin's numbering from the multi-layered cartesian
    # in anticipation. For PEC and dielectric sphere;
    # 1 = interior, 2 = exterior

    if r >= sp.radius
        return 2
    else
        return 1
    end
end

function wavenumber(sp::PECSphere, ex::Excitation, r)
    ε = ex.embedding.ε
    μ = ex.embedding.μ

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    if layer(sp, r) == 2
        return k
    else
        return typeof(k)(0.0)
    end
end

function wavenumber(sp::DielectricSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    return k
end

function impedance(sp::DielectricSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    return sqrt(μ / ε)
end

function impedance(sp::PECSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        return promote_type(typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
    end

    return sqrt(μ / ε)
end
