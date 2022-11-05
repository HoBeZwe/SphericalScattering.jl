
abstract type SphericalMode <: Excitation end

struct SphericalModeTE{T,R,C,In<:Integer} <: SphericalMode
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    m::In
    n::In
    c::In
    center::SVector{3,R}
    orientation::SVector{3,R}
end

struct SphericalModeTM{T,R,C,In<:Integer} <: SphericalMode
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    m::In
    n::In
    c::In
    center::SVector{3,R}
    orientation::SVector{3,R}
end

SphericalModeTE(;
    embedding   = Medium(ε0, μ0),
    wavenumber  = error("missing argument `wavenumber`"),
    amplitude   = 1.0,
    m           = 0,
    n           = 1,
    c           = 1,
    center      = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
    orientation = SVector{3,typeof(wavenumber)}(0.0, 0.0, 1.0),
) = SphericalModeTE(embedding, wavenumber, amplitude, m, n, c, center, orientation)

SphericalModeTM(;
    embedding   = Medium(ε0, μ0),
    wavenumber  = error("missing argument `wavenumber`"),
    amplitude   = 1.0,
    m           = 0,
    n           = 1,
    c           = 1,
    center      = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
    orientation = SVector{3,typeof(wavenumber)}(0.0, 0.0, 1.0),
) = SphericalModeTM(embedding, wavenumber, amplitude, m, n, c, center, orientation)




"""
    getFieldType(excitation::SphericalMode, quantity::ElectricField)

Do nothing.
"""
function getFieldType(excitation::SphericalMode, quantity::Field)
    return quantity, excitation
end



"""
    getFieldType(excitation::SphericalModeTE, quantity::MagneticField)

Exchange TE and TM. Other prefactor: im / ZF.
"""
function getFieldType(excitation::SphericalModeTE, quantity::MagneticField)

    exc = SphericalModeTM(
        excitation.embedding,
        excitation.wavenumber,
        im * sqrt(excitation.embedding.ε / excitation.embedding.μ) * excitation.amplitude,
        excitation.m,
        excitation.n,
        excitation.c,
        excitation.center,
        excitation.orientation,
    )
    return ElectricField(quantity.locations), exc
end



"""
    getFieldType(excitation::SphericalModeTM, quantity::MagneticField)

Exchange TE and TM. Other prefactor: im / ZF.
"""
function getFieldType(excitation::SphericalModeTM, quantity::MagneticField)

    exc = SphericalModeTE(
        excitation.embedding,
        excitation.wavenumber,
        im * sqrt(excitation.embedding.ε / excitation.embedding.μ) * excitation.amplitude,
        excitation.m,
        excitation.n,
        excitation.c,
        excitation.center,
        excitation.orientation,
    )
    return ElectricField(quantity.locations), exc
end
