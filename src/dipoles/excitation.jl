
abstract type Dipole <: Excitation end

struct HertzianDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    frequency::R
    amplitude::T
    center::SVector{3,R}
    orientation::SVector{3,R}

    # inner constructor: normalize orientation
    function HertzianDipole(
        embedding::Medium{C}, frequency::R, amplitude::T, center::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}
        or_normalized = normalize(orientation)
        new{T,R,C}(embedding, frequency, amplitude, center, or_normalized)
    end
end

struct FitzgeraldDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    frequency::R
    amplitude::T
    center::SVector{3,R}
    orientation::SVector{3,R}

    function FitzgeraldDipole(
        embedding::Medium{C}, frequency::R, amplitude::T, center::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}
        or_normalized = normalize(orientation)
        new{T,R,C}(embedding, frequency, amplitude, center, or_normalized)
    end
end

"""
    ex = HertzianDipole(
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            center      = SVector{3,typeof(frequency)}(0.0, 0.0, 0.0),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
HertzianDipole(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    center      = SVector{3,typeof(frequency)}(0.0, 0.0, 0.0),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = HertzianDipole(embedding, frequency, amplitude, center, orientation)


"""
    ex = FitzgeraldDipole(;
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            center      = SVector{3,typeof(frequency)}(0.0, 0.0, 0.0),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
FitzgeraldDipole(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    center      = SVector{3,typeof(frequency)}(0.0, 0.0, 0.0),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = FitzgeraldDipole(embedding, frequency, amplitude, center, orientation)






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
        exc = FitzgeraldDipole(embedding, excitation.frequency, -excitation.amplitude, excitation.center, excitation.orientation) # change sign (duality realations)
        return MagneticField(quantity.locations), exc

    elseif typeof(quantity) == MagneticField
        exc = FitzgeraldDipole(embedding, excitation.frequency, excitation.amplitude, excitation.center, excitation.orientation)
        return ElectricField(quantity.locations), exc

    else
        exc = FitzgeraldDipole(embedding, excitation.frequency, -excitation.amplitude, excitation.center, excitation.orientation) # change sign (duality realations)
        return quantity, exc
    end
end
