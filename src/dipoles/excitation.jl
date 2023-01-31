
abstract type Dipole <: Excitation end

struct HertzianDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
    orientation::SVector{3,R}

    # inner constructor: normalize orientation and check orientation
    function HertzianDipole(
        embedding::Medium{C}, frequency::R, amplitude::T, position::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}

        or_normalized = normalize(orientation)
        orientation × position == SVector{3,R}(0, 0, 0) || @info "The dipole is not placed perpendicular to a sphere."

        new{T,R,C}(embedding, frequency, amplitude, position, or_normalized)
    end
end

struct FitzgeraldDipole{T,R,C} <: Dipole
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
    orientation::SVector{3,R}

    # inner constructor: normalize orientation and check orientation
    function FitzgeraldDipole(
        embedding::Medium{C}, frequency::R, amplitude::T, position::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}

        or_normalized = normalize(orientation)
        orientation × position == SVector{3,R}(0, 0, 0) || @info "The dipole is not placed perpendicular to a sphere."

        new{T,R,C}(embedding, frequency, amplitude, position, or_normalized)
    end
end


"""
    ex = HertzianDipole(
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            position    = error("missing argument `position`"),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
HertzianDipole(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    position    = error("missing argument `position`"),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = HertzianDipole(embedding, frequency, amplitude, position, orientation)


"""
    ex = FitzgeraldDipole(;
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            position    = error("missing argument `position`"),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
FitzgeraldDipole(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    position    = error("missing argument `position`"),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = FitzgeraldDipole(embedding, frequency, amplitude, position, orientation)






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
        exc = FitzgeraldDipole(embedding, excitation.frequency, -excitation.amplitude, excitation.position, excitation.orientation) # change sign (duality realations)
        return MagneticField(quantity.locations), exc

    elseif typeof(quantity) == MagneticField
        exc = FitzgeraldDipole(embedding, excitation.frequency, excitation.amplitude, excitation.position, excitation.orientation)
        return ElectricField(quantity.locations), exc

    else
        exc = FitzgeraldDipole(embedding, excitation.frequency, -excitation.amplitude, excitation.position, excitation.orientation) # change sign (duality realations)
        return quantity, exc
    end
end
