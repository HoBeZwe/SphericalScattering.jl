
struct PlaneWave{T,R,C} <: Excitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    direction::SVector{3,R}
    polarization::SVector{3,R}

    # inner constructor: normalize direction and polarization and check their orthogonality
    function PlaneWave(
        embedding::Medium{C}, frequency::R, amplitude::T, direction::SVector{3,R}, polarization::SVector{3,R}
    ) where {T,R,C}

        dot(direction, polarization) == 0.0 ||
            error("No porper definition of a plane wave: direction and polarization vector are not perpendicular!")

        dir_normalized = normalize(direction)
        pol_normalized = normalize(polarization)

        new{T,R,C}(embedding, frequency, amplitude, dir_normalized, pol_normalized)
    end
end


"""
    ex = planeWave(;
            embedding    = Medium(ε0, μ0),
            frequency    = error("missing argument `frequency`"),
            amplitude    = 1.0,
            direction    = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
            polarization = SVector{3,typeof(frequency)}(1.0, 0.0, 0.0),
    )
"""
planeWave(;
    embedding    = Medium(ε0, μ0),
    frequency    = error("missing argument `frequency`"),
    amplitude    = 1.0,
    direction    = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
    polarization = SVector{3,typeof(frequency)}(1.0, 0.0, 0.0),
) = PlaneWave(embedding, frequency, amplitude, direction, polarization)

