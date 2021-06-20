
abstract type Dipole end

struct HertzianDipole{T} <: Dipole
    embedding::Medium{T}
    wavenumber::T
    amplitude::T
    center::SVector{3,T}
    orientation::SVector{3,T}
end

struct FitzgeraldDipole{T} <: Dipole
    embedding::Medium{T}
    wavenumber::T
    amplitude::T
    center::SVector{3,T}
    orientation::SVector{3,T}
end

HertzianDipole(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
center      = SVector(0.0,0.0,0.0),
orientation = SVector(0.0,0.0,1.0)
) = HertzianDipole(embedding, wavenumber, amplitude, center, orientation)

FitzgeraldDipole(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
center      = SVector(0.0,0.0,0.0),
orientation = SVector(0.0,0.0,1.0)
) = FitzgeraldDipole(embedding, wavenumber, amplitude, center, orientation)






"""
    getFieldType(excitation::HertzianDipole, quantity::Field)

Do nothing.
"""
function getFieldType(excitation::HertzianDipole, quantity::Field)
    return quantity, excitation
end


"""
    getFieldType(excitation::FitzgeraldDipole, quantity::Field)

Exchange electric and magnetic field for a magnetic current + apply duality relations.
"""
function getFieldType(excitation::FitzgeraldDipole, quantity::Field)
   
    embedding = Medium(excitation.embedding.μ, excitation.embedding.ε) # exchange μ and ε (duality relations)

    if typeof(quantity) == ElectricField
        exc = FitzgeraldDipole(embedding, excitation.wavenumber, -excitation.amplitude, excitation.center, excitation.orientation) # change sign (duality realations)
        return MagneticField(quantity.locations), exc

    elseif typeof(quantity) == MagneticField
        exc = FitzgeraldDipole(embedding, excitation.wavenumber, excitation.amplitude, excitation.center, excitation.orientation) 
        return ElectricField(quantity.locations), exc

    else
        exc = FitzgeraldDipole(embedding, excitation.wavenumber, -excitation.amplitude, excitation.center, excitation.orientation) # change sign (duality realations)
        return quantity, exc
    end
end