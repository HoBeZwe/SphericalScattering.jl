
"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

Compute the field scattered by a PEC sphere excited by a dipole at some position and orientation.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium
    excitation.orientation × excitation.position == SVector{3,T}(0, 0, 0) || error("The dipole is not perpendicular to the sphere.")

    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- rotate coordinates
    points = rotate(excitation, quantity.locations; inverse=true)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(points)
        F[ind] = scatteredfield(sphere, exc, point, fieldType; parameter=parameter)
    end

    # --- rotate resulting field
    rotate!(excitation, F; inverse=false)

    return F
end



"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = wavenumber(excitation)
    T  = typeof(k)
    Il = excitation.amplitude
    z0 = norm(excitation.position)   # distance Dipole-origin

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    eps = parameter.relativeAccuracy

    Er  = Complex{T}(0.0) # initialize
    Eϑ  = Complex{T}(0.0) # initialize
    δEr = T(Inf)
    n   = 0

    r = point_sph[1]

    kr  = k * r
    kz0 = k * z0
    ka  = k * sphere.radius

    r <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    try
        while δEr > eps || n < 10
            n += 1
            expansion_coeff = 2 * n + 1  # variable part of expansion coefficients

            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkr = hankelh2(n + T(0.5), kr)   # Hankel function 2nd kind (without sqrt-term)
            Hkz0 = hankelh2(n + T(0.5), kz0) # Hankel function 2nd kind (without sqrt-term)
            Hkrt = hankelh2(n + T(1.5), kr)  # for derivative of Riccati-Hankel function 2nd kind

            dHkr = (n + 1) * sqrt(z0 / r) * Jka_Hka * Hkr * Hkz0 - k * sqrt(r * z0) * Jka_Hka * Hkz0 * Hkrt  # HkR * derivative of Riccati-Hankel function 2nd kind

            ΔEr = expansion_coeff * n * (n + 1) * Jka_Hka * Hkr * Hkz0 * Pl(cosϑ, n)
            ΔEϑ = expansion_coeff * dHkr * dnPl(cosϑ, n, 1)

            isodd(n) && (δEr = abs(ΔEr / Er)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Er += ΔEr
            Eϑ += ΔEϑ
        end
    catch
        # print("did not converge: n=$n\n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    ZF = sqrt(μ / ε)

    Er *= +im * Il / 8 / k / z0 * sqrt(r / z0) / r / r * ZF
    Eϑ *= -im * Il / 8 / k / z0 / z0 * sinϑ / r * ZF

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, 0.0), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = wavenumber(excitation)
    T  = typeof(k)
    Il = excitation.amplitude
    z0 = norm(excitation.position)    # distance Dipole-origin

    eps = parameter.relativeAccuracy

    Hϕ = Complex{T}(0.0) # initialize
    δH = T(Inf)
    n = 0

    r = point_sph[1]

    kr  = k * r
    kz0 = k * z0
    ka  = k * sphere.radius

    r <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    try
        while δH > eps || n < 10
            n += 1
            expansion_coeff = 2 * n + 1 # variable part of expansion coefficients

            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + T(0.5), kz0) # Hankel function 2nd kind
            Hkr  = hankelh2(n + T(0.5), kr) # Hankel function 2nd kind

            ΔH = expansion_coeff * Jka_Hka * Hkz0 * Hkr * dnPl(cosϑ, n, 1)
            isodd(n) && (δH = abs(ΔH / Hϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Hϕ += ΔH
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Hϕ *= -Il / z0 / sqrt(z0 * r) * sinϑ / 8

    return SVector{3,Complex{T}}(-Hϕ * sin(point_sph[3]), Hϕ * cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::HertzianDipole, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::HertzianDipole, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = wavenumber(excitation)
    T  = typeof(k)
    Il = excitation.amplitude
    z0 = norm(excitation.position)    # distance dipole-origin

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    eps = parameter.relativeAccuracy

    Eϑ = Complex{T}(0.0) # initialize
    δE = T(Inf)
    n  = 0

    kz0 = k * z0
    ka  = k * sphere.radius

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    try
        while δE > eps || n < 10
            n += 1
            expansion_coeff = 2 * n + 1   # variable part of expansion coefficients

            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + T(0.5), kz0) # Hankel function 2nd kind

            ΔE = expansion_coeff * Jka_Hka * Hkz0 * im^(n + 1) * dnPl(cosϑ, n, 1)
            isodd(n) && (δE = abs(ΔE / Eϑ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϑ += ΔE
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϑ *= -Il / k * sinϑ / z0 / 4 * sqrt(k / z0 / 2 / π) * sqrt(μ / ε)

    return SVector{3,Complex{T}}(
        Eϑ * cos(point_sph[2]) * cos(point_sph[3]), Eϑ * cos(point_sph[2]) * sin(point_sph[3]), -Eϑ * sin(point_sph[2])
    ) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::FitzgeraldDipole, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::FitzgeraldDipole, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = wavenumber(excitation)
    T  = typeof(k)
    Ul = excitation.amplitude
    z0 = norm(excitation.position)    # distance dipole-origin

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    eps = parameter.relativeAccuracy

    Eϕ = Complex{T}(0.0) # initialize
    δE = T(Inf)
    n  = 0

    kz0 = k * z0
    ka  = k * sphere.radius

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    try
        while δE > eps || n < 10
            n += 1
            expansion_coeff = 2 * n + 1 # variable part of expansion coefficients

            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + T(0.5), kz0) # Hankel function 2nd kind

            ΔE = expansion_coeff * Jka_Hka * Hkz0 * im^(n + 1) * dnPl(cosϑ, n, 1)
            isodd(n) && (δE = abs(ΔE / Eϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϕ += ΔE
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϕ *= -Ul / k * sinϑ / z0 / 4 * sqrt(k / z0 / 2 / π)

    return SVector{3,Complex{T}}(-Eϕ * sin(point_sph[3]), Eϕ * cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatterCoeff(sphere::PECSphere, excitation::FitzgeraldDipole, n::Int, ka)

Compute scattering coefficient for Fitzgerald dipole.
"""
function scatterCoeff(sphere::PECSphere, excitation::FitzgeraldDipole, n::Int, ka)

    T = typeof(ka)
    Jka = besselj(n + T(0.5), ka)  # Bessel function 1st kind
    Hka = hankelh2(n + T(0.5), ka) # Hankel function 2nd kind

    return Jka / Hka
end



"""
    scatterCoeff(sphere::PECSphere, excitation::HertzianDipole, n::Int, ka)

Compute scattering coefficient for Hertzian dipole.
"""
function scatterCoeff(sphere::PECSphere, excitation::HertzianDipole, n::Int, ka)

    T = typeof(ka)
    Jka = (n + 1) * besselj(n + T(0.5), ka) - ka * besselj(n + T(1.5), ka)  # derivative spherical Bessel function 1st kind (without sqrt factor)
    Hka = (n + 1) * hankelh2(n + T(0.5), ka) - ka * hankelh2(n + T(1.5), ka) # derivative spherical Hankel function 2nd kind (without sqrt factor)

    return Jka / Hka
end
