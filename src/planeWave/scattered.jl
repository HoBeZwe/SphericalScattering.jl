
"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())
    
Compute the electric field scattered by a PEC sphere, for an incident plane wave.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    T = typeof(excitation.wavenumber)

    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- rotate coordinates
    # rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, point, quantity; parameter=parameter)
    end

    # --- rotate resulting field
    # rotate!(F, excitation.rotation)

    return F
end



"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, for an incident plane wave
travelling in +z-direction with E-field polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
    E₀ = excitation.amplitude

    T = typeof(k)

    eps = parameter.relativeAccuracy

    point_sph[1] <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0

    Er = Complex{T}(0.0)
    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    δE = T(Inf)
    n = 0

    kr = k * point_sph[1]
    ka = k * sphere.radius

    sinϑ = abs(sin(point_sph[2]))  # note: theta only defined from from 0 to pi
    cosϑ = cos(point_sph[2])       # ok for theta > pi
    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    # first two values of the Associated Legendre Polynomial
    plm = Vector{T}()
    push!(plm, -sinϑ)
    push!(plm, -T(3.0) * sinϑ * cosϑ)

    s = sqrt(π / 2 / kr)

    try
        while δE > eps
            n += 1

            aₙ, bₙ = scatterCoeff(sphere, excitation, n, ka)

            Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)

            ΔEr = +(cosϕ / (im * kr^2)) * aₙ * Nn_r
            ΔEϑ = -(cosϕ / kr) * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
            ΔEϕ = +(sinϕ / kr) * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)

            Er += ΔEr
            Eϑ += ΔEϑ
            Eϕ += ΔEϕ

            δE = (abs(ΔEr) + abs(ΔEϑ) + abs(ΔEϕ)) / (abs(Er) + abs(Eϑ) + abs(Eϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    return convertSpherical2Cartesian(E₀ .* SVector(Er, Eϑ, Eϕ), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field scattered by a PEC sphere, for an incident plane wave
travelling in +z-direction with E-field polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
    T = typeof(k)
    μ = excitation.embedding.μ
    ε = excitation.embedding.ε
    η = sqrt(μ / ε)
    H₀ = 1 / η * excitation.amplitude

    eps = parameter.relativeAccuracy

    point_sph[1] <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0

    Hr = Complex{T}(0.0)
    Hϑ = Complex{T}(0.0)
    Hϕ = Complex{T}(0.0)

    δH = T(Inf)
    n = 0

    kr = k * point_sph[1]
    ka = k * sphere.radius

    sinϑ = abs(sin(point_sph[2]))  # note: theta only defined from from 0 to pi
    cosϑ = cos(point_sph[2])       # ok for theta > pi
    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    # first two values of the Associated Legendre Polynomial
    plm = Vector{T}()
    push!(plm, -sinϑ)
    push!(plm, -T(3.0) * sinϑ * cosϑ)

    s = sqrt(π / 2 / kr)

    try
        while δH > eps
            n += 1

            aₙ, bₙ = scatterCoeff(sphere, excitation, n, ka)
            Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)

            ΔHr = +(sinϕ / (im * kr^2)) * bₙ * Nn_r
            ΔHϑ = -(sinϕ / kr) * (aₙ * Mn_ϑ + bₙ * Nn_ϑ)
            ΔHϕ = -(cosϕ / kr) * (aₙ * Mn_ϕ + bₙ * Nn_ϕ)

            Hr += ΔHr
            Hϑ += ΔHϑ
            Hϕ += ΔHϕ

            δH = (abs(ΔHr) + abs(ΔHϑ) + abs(ΔHϕ)) / (abs(Hr) + abs(Hϑ) + abs(Hϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    return convertSpherical2Cartesian(H₀ .* SVector(Hr, Hϑ, Hϕ), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the (electric) far-field scattered by a PEC sphere, for an incident plane wave
travelling in +z-direction with E-field polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
    E₀ = excitation.amplitude

    T = typeof(k)
    eps = parameter.relativeAccuracy

    sinϑ = abs(sin(point_sph[2]))  # note: theta only defined from from 0 to pi
    cosϑ = cos(point_sph[2])       # ok for theta > pi
    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    # first two values of the Associated Legendre Polynomial
    plm = Vector{T}()
    push!(plm, -sinϑ)
    push!(plm, -T(3.0) * sinϑ * cosϑ)

    ka = k * sphere.radius

    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    δE = T(Inf)
    n = 0

    try
        while δE > eps || n < 10
            n += 1

            aₙ, bₙ = scatterCoeff(sphere, excitation, n, ka)
            Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, ka, plm, cosϑ, sinϑ, n)

            # See Jin, (7.4.44)-(7.4.45)
            # We deviate from Jin by replacing exp(-im*kr) / kr with 1/k
            ΔEϑ = -im * cosϕ * 1 / k * im^n * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
            ΔEϕ = +im * sinϕ * 1 / k * im^n * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)

            Eϑ += ΔEϑ
            Eϕ += ΔEϕ

            δE = (abs(ΔEϑ) + abs(ΔEϕ)) / (abs(Eϑ) + abs(Eϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship associated Legendre polynomial
        end
    catch

    end

    return convertSpherical2Cartesian(E₀ .* SVector{3,Complex{T}}(0.0, Eϑ, Eϕ), point_sph)
end



"""
    scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int, ka)

Compute scattering coefficients for a plane wave travelling in -z direction with polarization in x-direction.
"""
function scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int, ka)

    T = typeof(excitation.wavenumber)
    s = sqrt(π / 2 / ka)

    Ĵ  = ka * s * besselj(n + T(0.5), ka)   # Riccati-Bessel function
    Ĥ  = ka * s * hankelh2(n + T(0.5), ka)  # Riccati-Hankel function
    Ĵ2 = ka * s * besselj(n - T(0.5), ka)   # for derivate needed
    Ĥ2 = ka * s * hankelh2(n - T(0.5), ka)  # for derivate needed

    # Use recurrence relationship
    dĴ = (Ĵ2 - n / ka * Ĵ)    # derivative Riccati-Bessel function
    dĤ = (Ĥ2 - n / ka * Ĥ)    # derivative Riccati-Hankel function

    aₙ = -im^(-T(n)) * (dĴ / dĤ) * (2 * n + 1) / (n * (n + 1))  # Jin (7.4.41)
    bₙ = -im^(-T(n)) * (Ĵ / Ĥ) * (2 * n + 1) / (n * (n + 1))    # Jin (7.4.42)

    return aₙ, bₙ
end



"""
    expansion(sphere::PECSphere, excitation::PlaneWave, plm, kr, s, cosϑ, p, n::Int) 

Compute functional dependencies of the Mie series for a plane wave travelling in -z direction with polarization in x-direction.
"""
function expansion(sphere::PECSphere, excitation::PlaneWave, plm, kr, s, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.wavenumber)

    Ĥ  = kr * s * hankelh2(n + T(0.5), kr)     # Riccati-Hankel function
    Ĥ2 = kr * s * hankelh2(n - T(0.5), kr)

    dĤ = Ĥ2 - n / kr * Ĥ          # derivative of Riccati-Hankel function

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = Ĥ * p / sinϑ
        Mn_ϕ = Ĥ * dp

        Nn_ϑ = dĤ * dp
        Nn_ϕ = dĤ * p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -Ĥ * aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -dĤ * aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = Ĥ * aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -dĤ * aux
        Nn_ϕ = -Nn_ϑ
    end

    Nn_r = n * (n + 1) * Ĥ * p
    Nn_ϑ = im * Nn_ϑ
    Nn_ϕ = im * Nn_ϕ

    return Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end



"""
    expansion(sphere::PECSphere, excitation::PlaneWave, ka, plm, cosϑ, sinϑ, n::Int)

Compute far-field functional dependencies of the Mie series for a plane wave travelling in -z direction with polarization in x-direction.
"""
function expansion(sphere::PECSphere, excitation::PlaneWave, ka, plm, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.wavenumber)

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = p / sinϑ
        Mn_ϕ = dp

        Nn_ϑ = dp
        Nn_ϕ = p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -aux
        Nn_ϕ = -Nn_ϑ
    end

    return Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end
