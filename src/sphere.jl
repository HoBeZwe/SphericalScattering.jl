
struct Medium{T}
    ε::T
    μ::T
end

abstract type Sphere end

struct DielectricSphere{T} <: Sphere
    radius::T
    embedding::Medium{T}
    filling::Medium{T}
end

struct PECSphere{T} <: Sphere
    radius::T
    embedding::Medium{T}
end

PECSphere(;
radius      = error("missing argument `wavenumber`"),
embedding   = Medium(ε0, μ0)
) = PECSphere(radius, embedding)

struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{N,Medium{C}}
    embedding::Medium{C}
end

LayeredSphere(;
radii = error("Missing argument `radii`"),
filling = error("`missing argument `filling`"),
embedding = Medium(ε0, μ0)
) = LayeredSphere(radii, filling, embedding)

struct LayeredSpherePEC{N,D,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{D,Medium{C}}
    embedding::Medium{C}
end

LayeredSpherePEC(;
    radii = error("Missing argument `radii`"),
    filling = error("Missing argument `filling`"),
    embedding = Medium(ε0, μ0)
) = LayeredSpherePEC(radii, filling, embedding)
