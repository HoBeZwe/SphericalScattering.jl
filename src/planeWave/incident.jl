
"""
    field(excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field of a plane wave.
"""
function field(excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

    F = zeros(SVector{3,Complex{Float64}}, length(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity, parameter=parameter)
    end

    return F
end



"""
    field(excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field of a plane wave.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = excitation.wavenumber

    d = excitation.direction
    p = excitation.polarization

    return a * exp(-im * k * dot(d, point)) * p
end



"""
    field(excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the magnetic field of a plane wave.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = excitation.wavenumber

    d = excitation.direction
    p = excitation.polarization

    ε  = excitation.embedding.ε
    μ  = excitation.embedding.μ

    return a * sqrt(ε / μ) * exp(-im * k * dot(d, point)) * (d × p)
end