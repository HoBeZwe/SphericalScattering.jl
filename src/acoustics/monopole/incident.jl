
"""
    field(excitation::AcousticMonopole, quantity::Field; parameter::Parameter=Parameter())

Compute the field of an acoustic monopole.
"""
function field(excitation::AcousticMonopole, quantity::Pressure; parameter::Parameter=Parameter(), zeroRadius=0.0)

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
    field(excitation::AcousticMonopole, point, quantity::Pressure; parameter::Parameter=Parameter())

Compute the field of an acoustic monopole.

The point is in Cartesian coordinates.
"""
function field(excitation::AcousticMonopole, point, quantity::Pressure; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    r = norm(point - excitation.position)

    return a * cis(-k * r)
end