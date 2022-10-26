
struct Medium{C}
    ε::C
    μ::C
end

abstract type Sphere end

struct DielectricSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
    filling::Medium{C}
end

DielectricSphere(;
radius      = error("missing argument `wavenumber`"),
embedding   = Medium(ε0, μ0),
filling = error("missing argument `filling`")
) = DielectricSphere(radius, embedding, filling)

struct PECSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
end

PECSphere(;
radius      = error("missing argument `wavenumber`"),
embedding   = Medium(ε0, μ0)
) = PECSphere(radius, embedding)

struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    embedding::Medium{C}
    filling::SVector{N,Medium{C}}
end

LayeredSphere(;
radii = error("Missing argument `radii`"),
embedding = Medium(ε0, μ0),
filling = error("`missing argument `filling`")
) = LayeredSphere(radii, embedding, filling)

struct LayeredSpherePEC{N,D,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{D,Medium{C}}
    embedding::Medium{C}
end

LayeredSpherePEC(;
    radii = error("Missing argument `radii`"),
    embedding = Medium(ε0, μ0),
    filling = error("Missing argument `filling`")
) = LayeredSpherePEC(radii, embedding, filling)
