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




abstract type Sphere end

struct DielectricSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
    filling::Medium{C}
end

function DielectricSphere(r::R, embedding::Medium{C1}, filling::Medium{C2}) where {R,C1,C2}

    C = promote_type(C1, C2)

    DielectricSphere(r, Medium{C}(embedding), Medium{C}(filling))
end


"""
    DielectricSphere(
        radius      = error("missing argument `radius`"), 
        embedding   = Medium(ε0, μ0), 
        filling     = error("missing argument `filling`")
    )

Constructor for the dielectric sphere.
"""
DielectricSphere(; radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0), filling=error("missing argument `filling`")) =
    DielectricSphere(radius, embedding, filling)



struct DielectricSphereThinImpedanceLayer{R,C} <: Sphere
    radius::R
    thickness::R
    embedding::Medium{C}
    thinlayer::Medium{C}
    filling::Medium{C}
end

function DielectricSphereThinImpedanceLayer(
    r::R1, d::R2, embedding::Medium{C1}, thinlayer::Medium{C2}, filling::Medium{C3}
) where {R1,R2,C1,C2,C3}

    R = promote_type(R1, R2)
    C = promote_type(C1, C2, C3)

    DielectricSphereThinImpedanceLayer(R(r), R(d), Medium{C}(embedding), Medium{C}(thinlayer), Medium{C}(filling))
end

"""
    DielectricSphereThinImpedanceLayer(
        radius      = error("missing argument `radius`"),
        thickness   = error("missing argument `thickness` of the coating"),
        embedding   = Medium(ε0, μ0),
        thinlayer   = error("missing argument `thinlayer`"),
        filling     = error("missing argument `filling`")
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
    embedding=Medium(ε0, μ0),
    thinlayer=error("missing argument `thinlayer`"),
    filling=error("missing argument `filling`"),
) = DielectricSphereThinImpedanceLayer(radius, thickness, embedding, thinlayer, filling)

struct PECSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
end

"""
    PECSphere( 
        radius      = error("missing argument `radius`"), 
        embedding   = Medium(ε0, μ0)
    )

Constructor for the PEC sphere.
"""
PECSphere(; radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0)) = PECSphere(radius, embedding)



struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    embedding::Medium{C}
    filling::SVector{N,Medium{C}}
end

"""
    LayeredSphere( 
        radii       = error("Missing argument `radii`"), 
        embedding   = Medium(ε0, μ0), 
        filling     = error("`missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
LayeredSphere(; radii=error("Missing argument `radii`"), embedding=Medium(ε0, μ0), filling=error("`missing argument `filling`")) =
    LayeredSphere(radii, embedding, filling)




struct LayeredSpherePEC{N,D,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{D,Medium{C}}
    embedding::Medium{C}
end

"""
    LayeredSpherePEC( 
        radii       = error("Missing argument `radii`"), 
        embedding   = Medium(ε0, μ0), 
        filling     = error("Missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
LayeredSpherePEC(; radii=error("Missing argument `radii`"), embedding=Medium(ε0, μ0), filling=error("Missing argument `filling`")) =
    LayeredSpherePEC(radii, filling, embedding)




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
    ε = sp.embedding.ε
    μ = sp.embedding.μ

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
        ε = sp.embedding.ε
        μ = sp.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    return k
end

function impedance(sp::DielectricSphere, r)
    if layer(sp, r) == 2
        ε = sp.embedding.ε
        μ = sp.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    return sqrt(μ / ε)
end

function impedance(sp::PECSphere, r)
    if layer(sp, r) == 2
        ε = sp.embedding.ε
        μ = sp.embedding.μ
    else
        return promote_type(typeof(sp.embedding.ε), typeof(sp.embedding.μ))(0.0)
    end

    return sqrt(μ / ε)
end
