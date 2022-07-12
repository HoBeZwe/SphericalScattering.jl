
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

struct PECSphere{C,R} <: Sphere
    radius::R
    embedding::Medium{C}
end

PECSphere(;
radius      = error("missing argument `wavenumber`"),
embedding   = Medium(ε0, μ0)
) = PECSphere(radius, embedding)