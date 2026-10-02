


struct AcousticPlaneWave{T,R,C} <: Excitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    direction::SVector{3,R}

    # inner constructor: normalize direction
    function AcousticPlaneWave(embedding::Medium{C}, frequency::R, amplitude::T, direction::SVector{3,R}) where {T,R,C}

        dir_normalized = normalize(direction)

        new{T,R,C}(embedding, frequency, amplitude, dir_normalized)
    end
end
