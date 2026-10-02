
"""
    Acoustic

Submodule bundling the constructors of the acoustic excitations, so that, e.g., `Acoustic.planeWave`
does not clash with the electromagnetic `planeWave`.
"""
module Acoustic

using ..SphericalScattering
using StaticArrays

"""
    ex = planeWave(;
            embedding = Medium(ε0, μ0),
            frequency = error("missing argument `frequency`"),
            amplitude = 1.0,
            direction = SVector{3,typeof(frequency)}(0.0, 0.0, 1.0)
    )

Constructor for an acoustic plane wave travelling in `direction`.
"""
planeWave(;
    embedding=SphericalScattering.Medium(ε0, μ0),
    frequency=error("missing argument `frequency`"),
    amplitude=1.0,
    direction=SVector{3,typeof(frequency)}(0.0, 0.0, 1.0),
) = SphericalScattering.AcousticPlaneWave(embedding, frequency, amplitude, direction)


"""
    ex = monopole(;
            position  = error("missing argument `position`"),
            frequency = error("missing argument `frequency`"),
            embedding = Medium(ε0, μ0),
            amplitude = 1.0
    )

Constructor for an acoustic monopole (point source) located at `position`.
"""
monopole(;
    position=error("missing argument `position`"),
    frequency=error("missing argument `frequency`"),
    embedding=SphericalScattering.Medium(ε0, μ0),
    amplitude=1.0,
) = SphericalScattering.AcousticMonopole(embedding, frequency, amplitude, position)

end
