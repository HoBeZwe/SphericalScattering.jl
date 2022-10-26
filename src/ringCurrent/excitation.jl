
abstract type RingCurrent <: Excitation end

struct ElectricRingCurrent{T,R,C} <: RingCurrent
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    radius::R
    center::SVector{3,R}
    rotation::SVector{2,R}
end

struct MagneticRingCurrent{T,R,C} <: RingCurrent
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    radius::R
    center::SVector{3,R}
    rotation::SVector{2,R}
end

electricRingCurrent(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
radius      = error("missing argument `radius`"),
center      = SVector{3,typeof(wavenumber)}(0.0,0.0,0.0),
rotation    = SVector{2,typeof(wavenumber)}(0.0,0.0)
) = ElectricRingCurrent(embedding, wavenumber, amplitude, radius, center, rotation)

magneticRingCurrent(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
radius      = error("missing argument `radius`"),
center      = SVector{3,typeof(wavenumber)}(0.0,0.0,0.0),
rotation    = SVector{2,typeof(wavenumber)}(0.0,0.0)
) = MagneticRingCurrent(embedding, wavenumber, amplitude, radius, center, rotation)






"""
    getFieldType(excitation::ElectricRingCurrent, quantity::Field)

Do nothing.
"""
function getFieldType(excitation::ElectricRingCurrent, quantity::Field)
    return quantity, excitation
end


"""
    getFieldType(excitation::MagneticRingCurrent, quantity::Field)

Exchange electric and magnetic field for a magnetic current.
"""
function getFieldType(excitation::MagneticRingCurrent, quantity::Field)
   
    embedding = Medium(excitation.embedding.μ, excitation.embedding.ε) # exchange μ and ε (duality relations)

    if typeof(quantity) == ElectricField
        exc = MagneticRingCurrent(embedding, excitation.wavenumber, -excitation.amplitude, excitation.radius, excitation.center, excitation.rotation) # change sign (duality realations)
        return MagneticField(quantity.locations), exc

    elseif typeof(quantity) == MagneticField
        exc = MagneticRingCurrent(embedding, excitation.wavenumber, excitation.amplitude, excitation.radius, excitation.center, excitation.rotation) 
        return ElectricField(quantity.locations), exc

    else
        exc = MagneticRingCurrent(embedding, excitation.wavenumber, -excitation.amplitude, excitation.radius, excitation.center, excitation.rotation) # change sign (duality realations)
        return quantity, exc
    end
end

