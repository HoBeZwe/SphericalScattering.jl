
abstract type SphericalMode <: Excitation end 

struct SphericalModeTE{T,C,In<:Integer} <: SphericalMode
    embedding::Medium{T}
    wavenumber::T
    amplitude::C
    m::In
    n::In
    s::In
    center::SVector{3,T}
    orientation::SVector{3,T}
end

struct SphericalModeTM{T,C,In<:Integer} <: SphericalMode
    embedding::Medium{T}
    wavenumber::T
    amplitude::C
    m::In
    n::In
    s::In
    center::SVector{3,T}
    orientation::SVector{3,T}
end

SphericalModeTE(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
m           = 0,
n           = 1,
s           = 1,
center      = SVector(0.0,0.0,0.0),
orientation = SVector(0.0,0.0,1.0)
) = SphericalModeTE(embedding, wavenumber, amplitude, m, n, s, center, orientation)

SphericalModeTM(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
m           = 0,
n           = 1,
s           = 1,
center      = SVector(0.0,0.0,0.0),
orientation = SVector(0.0,0.0,1.0)
) = SphericalModeTM(embedding, wavenumber, amplitude, m, n, s, center, orientation)