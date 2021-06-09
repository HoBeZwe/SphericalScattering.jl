
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