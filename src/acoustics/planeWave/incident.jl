
"""
    field(excitation::AcousticPlaneWave, quantity::Field; parameter::Parameter=Parameter())

Compute the field of an acoustic plane wave.
"""
function field(excitation::AcousticPlaneWave, quantity::Pressure; parameter::Parameter=Parameter(), zeroRadius=0.0)

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        norm(point) < zeroRadius && continue
        F[ind] = field(excitation, point, quantity; parameter=parameter)
    end

    return F
end



"""
    field(excitation::AcousticPlaneWave, point, quantity::Pressure; parameter::Parameter=Parameter())

Compute the field of an acoustic plane wave.

The point is in Cartesian coordinates.
"""
function field(excitation::AcousticPlaneWave, point, quantity::Pressure; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    d = excitation.direction

    return a * cis(-k * dot(d, point))
end



"""
    field(excitation::AcousticPlaneWave, quantity::FarField; parameter::Parameter=Parameter())

Throw error since the far-field of a plane wave is not defined.
"""
function field(excitation::AcousticPlaneWave, quantity::FarField; parameter::Parameter=Parameter(), zeroRadius=0.0)

    return error("The far-field of a plane wave is not defined.")
end
