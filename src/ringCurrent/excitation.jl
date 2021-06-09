
struct ElectricRingCurrent{T}
    embedding::Medium{T}
    wavenumber::T
    amplitude::T
    radius::T
    center::SVector{3,T}
    rotation::SVector{2,T}
end

struct MagneticRingCurrent{T}
    embedding::Medium{T}
    wavenumber::T
    amplitude::T
    radius::T
    center::SVector{3,T}
    rotation::SVector{2,T}
end

electricRingCurrent(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
radius      = error("missing argument `radius`"),
center      = SVector(0.0,0.0,0.0),
rotation    = SVector(0.0,0.0)
) = ElectricRingCurrent(embedding, wavenumber, amplitude, radius, center, rotation)