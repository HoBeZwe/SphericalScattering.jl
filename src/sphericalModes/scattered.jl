
"""
    scatteredfield(sphere::PECSphere, excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a dipole at some position and orientation.
"""
function scatteredfield(sphere::PECSphere, excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    T = typeof(excitation.frequency)

    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- translate/rotate coordinates
    points = quantity.locations # translate(quantity.locations, -excitation.center)
    # rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(points)
        F[ind] = scatteredfield(sphere, exc, point, fieldType; parameter=parameter)
    end

    # --- rotate resulting field
    # rotate!(F, excitation.rotation)

    return F
end



"""
    scatteredfield(sphere::PECSphere, excitation::SphericalModeTE, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, when excited by a spherical mode travelling towards the origin.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::SphericalMode, point, quantity::Field; parameter::Parameter=Parameter())

    k  = wavenumber(excitation)
    ka = k * sphere.radius

    T = typeof(excitation)
    Q = typeof(quantity)

    # outward travelling wave
    Escat = T(
        excitation.embedding,
        excitation.frequency,
        excitation.amplitude,
        excitation.m,
        excitation.n,
        2,
        excitation.center,
        excitation.orientation,
    )

    # scatter coefficient
    γ = scatterCoeff(sphere, excitation, excitation.n, ka)
    E = field(Escat, Q([point]))

    return γ * E[1]
end



"""
    scatterCoeff(sphere::PECSphere, excitation::SphericalModeTE, n::Int, ka)

Compute scattering coefficients for a spherical TE mode travelling towards the origin.
"""
function scatterCoeff(sphere::PECSphere, excitation::SphericalModeTE, n::Int, ka)
    T = typeof(ka)
    return -hankelh2(n + T(0.5), ka) / hankelh1(n + T(0.5), ka)
end



"""
    scatterCoeff(sphere::PECSphere, excitation::SphericalModeTM, n::Int, ka)

Compute scattering coefficients for a spherical TM mode travelling towards the origin.
"""
function scatterCoeff(sphere::PECSphere, excitation::SphericalModeTM, n::Int, ka)

    T = typeof(ka)

    dH1 = (n + 1) * hankelh1(n + T(0.5), ka) - ka * hankelh1(n + T(1.5), ka)
    dH2 = (n + 1) * hankelh2(n + T(0.5), ka) - ka * hankelh2(n + T(1.5), ka)

    #dH1 = (n + 1) * besselj(n+0.5, ka) - ka * besselj(n+1.5, ka)
    #dH2 = (n + 1) * hankelh2(n+0.5, ka) - ka * hankelh2(n+1.5, ka)

    return -dH2 / dH1
end
