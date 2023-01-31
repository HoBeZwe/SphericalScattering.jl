
abstract type RingCurrent <: Excitation end

struct ElectricRingCurrent{T,R,C} <: RingCurrent
    embedding::Medium{C}
    frequency::R
    amplitude::T
    radius::R
    center::SVector{3,R}
    orientation::SVector{3,R}

    # inner constructor: normalize orientation and check orientation
    function ElectricRingCurrent(
        embedding::Medium{C}, frequency::R, amplitude::T, radius::R, center::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}

        or_normalized = normalize(orientation)
        orientation × center == SVector{3,R}(0, 0, 0) || @info "The ring current is not placed perpendicular to a sphere."

        new{T,R,C}(embedding, frequency, amplitude, radius, center, or_normalized)
    end
end

struct MagneticRingCurrent{T,R,C} <: RingCurrent
    embedding::Medium{C}
    frequency::R
    amplitude::T
    radius::R
    center::SVector{3,R}
    orientation::SVector{3,R}

    # inner constructor: normalize orientation and check orientation
    function MagneticRingCurrent(
        embedding::Medium{C}, frequency::R, amplitude::T, radius::R, center::SVector{3,R}, orientation::SVector{3,R}
    ) where {T,R,C}

        or_normalized = normalize(orientation)
        orientation × center == SVector{3,R}(0, 0, 0) || @info "The ring current is not placed perpendicular to a sphere."

        new{T,R,C}(embedding, frequency, amplitude, radius, center, or_normalized)
    end
end


"""
    ex = electricRingCurrent(;
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            radius      = error("missing argument `radius`"),
            center      = error("missing argument `center`"),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
electricRingCurrent(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    radius      = error("missing argument `radius`"),
    center      = error("missing argument `center`"),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = ElectricRingCurrent(embedding, frequency, amplitude, radius, center, orientation)


"""
    ex = magneticRingCurrent(;
            embedding   = Medium(ε0, μ0),
            frequency   = error("missing argument `frequency`"),
            amplitude   = 1.0,
            radius      = error("missing argument `radius`"),
            center      = error("missing argument `center`"),
            orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    )
"""
magneticRingCurrent(;
    embedding   = Medium(ε0, μ0),
    frequency   = error("missing argument `frequency`"),
    amplitude   = 1.0,
    radius      = error("missing argument `radius`"),
    center      = error("missing argument `center`"),
    orientation = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = MagneticRingCurrent(embedding, frequency, amplitude, radius, center, orientation)






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
        exc = MagneticRingCurrent(
            embedding, excitation.frequency, -excitation.amplitude, excitation.radius, excitation.center, excitation.orientation
        ) # change sign (duality realations)
        return MagneticField(quantity.locations), exc

    elseif typeof(quantity) == MagneticField
        exc = MagneticRingCurrent(
            embedding, excitation.frequency, excitation.amplitude, excitation.radius, excitation.center, excitation.orientation
        )
        return ElectricField(quantity.locations), exc

    else
        exc = MagneticRingCurrent(
            embedding, excitation.frequency, -excitation.amplitude, excitation.radius, excitation.center, excitation.orientation
        ) # change sign (duality realations)
        return quantity, exc
    end
end

