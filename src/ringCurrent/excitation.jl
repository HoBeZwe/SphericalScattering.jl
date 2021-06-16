
abstract type RingCurrent end

struct ElectricRingCurrent{T} <: RingCurrent
    embedding::Medium{T}
    wavenumber::T
    amplitude::T
    radius::T
    center::SVector{3,T}
    rotation::SVector{2,T}
end

struct MagneticRingCurrent{T} <: RingCurrent
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

magneticRingCurrent(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
radius      = error("missing argument `radius`"),
center      = SVector(0.0,0.0,0.0),
rotation    = SVector(0.0,0.0)
) = MagneticRingCurrent(embedding, wavenumber, amplitude, radius, center, rotation)






"""
    getFieldType(excitation::ElectricRingCurrent, quantity::Field)

Do nothing.
"""
function getFieldType(excitation::ElectricRingCurrent, quantity::Field)
    return quantity
end


"""
    getFieldType(excitation::MagneticRingCurrent, quantity::Field)

Exchange electric and magnetic field for a magnetic current.
"""
function getFieldType(excitation::MagneticRingCurrent, quantity::Field)
   
    embedding = Medium(excitation.embedding.μ, excitation.embedding.ε) # exchange μ and ε (duality relations)

    if typeof(quantity) == ElectricField
        excitation = MagneticRingCurrent(embedding, excitation.wavenumber, excitation.amplitude, excitation.radius, excitation.center, excitation.rotation)
        return MagneticField(quantity.locations)
    else
        excitation = MagneticRingCurrent(embedding, excitation.wavenumber, -excitation.amplitude, excitation.radius, excitation.center, excitation.rotation) # change sign (duality realations)
        return ElectricField(quantity.locations)
    end
end

