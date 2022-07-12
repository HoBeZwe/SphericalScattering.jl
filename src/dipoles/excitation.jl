
abstract type Dipole <: Excitation end 

struct HertzianDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    center::SVector{3,R}
    orientation::SVector{3,R}
end

struct FitzgeraldDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    center::SVector{3,R}
    orientation::SVector{3,R}
end

HertzianDipole(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
center      = SVector{3,typeof(wavenumber)}(0.0,0.0,0.0),
orientation = SVector{3,typeof(wavenumber)}(0.0,0.0,1.0)
) = HertzianDipole(embedding, wavenumber, amplitude, center, orientation)

FitzgeraldDipole(;
embedding   = Medium(ε0, μ0),
wavenumber  = error("missing argument `wavenumber`"),
amplitude   = 1.0,
center      = SVector{3,typeof(wavenumber)}(0.0,0.0,0.0),
orientation = SVector{3,typeof(wavenumber)}(0.0,0.0,1.0)
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