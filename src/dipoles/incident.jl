
"""
    field(excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field radiated by a magnetic/electric ring current at some position and with some orientation.
"""
function field(excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

    F = zeros(SVector{3,Complex{Float64}}, length(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- translate/rotate coordinates
    #points = translate(quantity.locations, excitation.center)
    # rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(exc, point, fieldType, parameter=parameter)
    end

    # --- rotate resulting field
    # rotate!(F, excitation.rotation)

    return F
end



"""
    field(excitation::HertzianDipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field radiated by a Hertzian dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = excitation.wavenumber
    ε  = excitation.embedding.ε
    μ  = excitation.embedding.μ

    r0 = excitation.center
    p  = excitation.orientation

    d = point - r0
    r = norm(d)
    n = d / r

    # r < 0.9 && return SVector(NaN, NaN, NaN)

    return Il / (4 * π) * sqrt(μ / ε) * exp(-im * k * r) * (k / r * cross(cross(n, p), n) + (1 / k / r^3 + im / r^2) * (3 * n * dot(n, p) - p))
end



"""
    field(excitation::HertzianDipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field radiated by a Hertzian dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = excitation.wavenumber

    r0 = excitation.center
    p  = excitation.orientation

    d = point - r0
    r = norm(d)
    n = d / r

    return -im * Il / (4 * π) * cross(n, p) * exp(-im * k * r) / r * (k + 1 / (im * r))
end
