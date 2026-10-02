


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


module Acoustic

using ..SphericalScattering
using StaticArrays

"""
    ex = planeWave(;
            embedding    = Medium(ε0, μ0),
            frequency    = error("missing argument `frequency`"),
            amplitude    = 1.0,
            direction    = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0)
    )
"""
planeWave(;
    embedding=SphericalScattering.Medium(ε0, μ0),
    frequency=error("missing argument `frequency`"),
    amplitude=1.0,
    direction=SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = SphericalScattering.AcousticPlaneWave(embedding, frequency, amplitude, direction)

end
