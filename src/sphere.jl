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

function wavenumber(sp::Sphere, ex::Excitation, r)
    if r >= sp.radius
        ε = sp.embedding.ε
        μ = sp.embedding.μ
    else
        if sp isa DielectricSphere
            ε = sp.filling.ε
            μ = sp.filling.μ
        elseif sp isa PECSphere
            return typeof(ex.frequency)(0.0)
        end
    end

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    return k
end

function impedance(sp::Sphere, r)
    if r >= sp.radius
        ε = sp.embedding.ε
        μ = sp.embedding.μ
    else
        if sp isa DielectricSphere
            ε = sp.filling.ε
            μ = sp.filling.μ
        elseif sp isa PECSphere
            return typeof(ex.frequency)(0.0)
        end
    end

    return sqrt(μ / ε)
end

"""
    DielectricSphere(;
        radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0), filling=error("missing argument `filling`")
    )

Constructor for the dielectric sphere.
"""
DielectricSphere(; radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0), filling=error("missing argument `filling`")) =
    DielectricSphere(radius, embedding, filling)


struct PECSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
end

"""
    PECSphere(; radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0))

Constructor for the PEC sphere.
"""
PECSphere(; radius=error("missing argument `radius`"), embedding=Medium(ε0, μ0)) = PECSphere(radius, embedding)


struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    embedding::Medium{C}
    filling::SVector{N,Medium{C}}
end

"""
    LayeredSphere(; radii=error("Missing argument `radii`"), embedding=Medium(ε0, μ0), filling=error("`missing argument `filling`"))

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
    LayeredSpherePEC(; radii=error("Missing argument `radii`"), embedding=Medium(ε0, μ0), filling=error("Missing argument `filling`"))

Constructor for the layered dielectric sphere.
"""
LayeredSpherePEC(; radii=error("Missing argument `radii`"), embedding=Medium(ε0, μ0), filling=error("Missing argument `filling`")) =
    LayeredSpherePEC(radii, filling, embedding)
