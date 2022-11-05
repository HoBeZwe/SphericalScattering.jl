
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

Compute the electric field scattered by a PEC sphere, for an incident plane wave travelling in -z direction with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
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

            An, Bn = scatterCoeff(sphere, excitation, n, ka)
            Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)

            ΔEr = Bn * Nn_r * cosϕ
            ΔEϑ = An * Mn_ϑ * cosϕ + Bn * Nn_ϑ * cosϕ
            ΔEϕ = -An * Mn_ϕ * sinϕ - Bn * Nn_ϕ * sinϕ

            Er += ΔEr
            Eϑ += ΔEϑ
            Eϕ += ΔEϕ

            δE = (abs(ΔEr) + abs(ΔEϑ) + abs(ΔEϕ)) / (abs(Er) + abs(Eϑ) + abs(Eϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    #return SVector(Er, Eϑ, Eϕ)
    return convertSpherical2Cartesian(SVector(Er, Eϑ, Eϕ), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field scattered by a PEC sphere, for an incident plane wave travelling in -z direction with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
    T = typeof(k)
    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

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

            An, Bn = scatterCoeff(sphere, excitation, n, ka)
            Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)

            ΔHr = An * Nn_r * sinϕ
            ΔHϑ = -Bn * Mn_ϑ * sinϕ + An * Nn_ϑ * sinϕ
            ΔHϕ = -Bn * Mn_ϕ * cosϕ + An * Nn_ϕ * cosϕ

            Hr += ΔHr
            Hϑ += ΔHϑ
            Hϕ += ΔHϕ

            δH = (abs(ΔHr) + abs(ΔHϑ) + abs(ΔHϕ)) / (abs(Hr) + abs(Hϑ) + abs(Hϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    return convertSpherical2Cartesian((im * sqrt(ε / μ)) .* SVector(Hr, Hϑ, Hϕ), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the (electric) far-field scattered by a PEC sphere, for an incident plane wave travelling in -z direction with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::PlaneWave, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k = excitation.wavenumber
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

    Sϑ = Complex{T}(0.0)
    Sϕ = Complex{T}(0.0)

    δS = T(Inf)
    n = 0

    try
        while δS > eps || n < 10
            n += 1

            ΔSϑ, ΔSϕ = expansion(sphere, excitation, ka, plm, cosϑ, sinϑ, n)

            Sϑ += ΔSϑ
            Sϕ += ΔSϕ

            δS = (abs(ΔSϑ) + abs(ΔSϕ)) / (abs(Sϑ) + abs(Sϑ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship associated Legendre polynomial
        end
    catch

    end

    Eϑ = Sϑ * cosϕ / k # Ruck, et. al. (3.1-5)
    Eϕ = -Sϕ * sinϕ / k

    #return SVector(0.0, Eϑ, Eϕ)
    return convertSpherical2Cartesian(SVector{3,Complex{T}}(0.0, Eϑ, Eϕ), point_sph)
end



"""
    scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int, ka)

Compute scattering coefficients for a plane wave travelling in -z direction with polarization in x-direction.
"""
function scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int, ka)

    T = typeof(excitation.wavenumber)
    s = sqrt(π / 2 / ka)

    J  = s * besselj(n + T(0.5), ka)   # spherical Bessel function
    H  = s * hankelh2(n + T(0.5), ka)  # spherical Hankel function
    J2 = s * besselj(n - T(0.5), ka)
    H2 = s * hankelh2(n - T(0.5), ka)

    # [k₀ * a * j_n(k₀ a)]'
    kaJ1P = (ka * J2 - n * J)    # derivatives spherical Bessel functions
    # [k₀ * a * h2_n(k₀ a)]'
    kaH1P = (ka * H2 - n * H)    # derivatives spherical Hankel functions

    An = -((im)^n) * (J / H) * (2 * n + 1) / (n * (n + 1))                # Ruck, et. al. (3.2-1)
    Bn = ((im)^(n + 1)) * (kaJ1P / kaH1P) * (2 * n + 1) / (n * (n + 1))     # Ruck, et. al. (3.2-2)

    return An, Bn
end



"""
    expansion(sphere::PECSphere, excitation::PlaneWave, plm, kr, s, cosϑ, p, n::Int) 

Compute functional dependencies of the Mie series for a plane wave travelling in -z direction with polarization in x-direction.
"""
function expansion(sphere::PECSphere, excitation::PlaneWave, plm, kr, s, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.wavenumber)

    Hn  = s * hankelh2(n + T(0.5), kr)     # spherical Hankel functions
    H2n = s * hankelh2(n - T(0.5), kr)

    krH1Pn = kr * H2n - n * Hn          # derivatives spherical Hankel functions

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = Hn * p / sinϑ
        Mn_ϕ = Hn * dp

        Nn_ϑ = krH1Pn * dp / kr
        Nn_ϕ = krH1Pn * p / (kr * sinϑ)

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -Hn * aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -krH1Pn * aux / kr
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = Hn * aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -krH1Pn * aux / kr
        Nn_ϕ = -Nn_ϑ
    end

    Nn_r = (n * (n + 1) / kr) * Hn * p

    return Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end



"""
    expansion(sphere::PECSphere, excitation::PlaneWave, ka, plm, cosϑ, sinϑ, n::Int)

Compute far-field functional dependencies of the Mie series for a plane wave travelling in -z direction with polarization in x-direction.
"""
function expansion(sphere::PECSphere, excitation::PlaneWave, ka, plm, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.wavenumber)
    An, Bn = scatterCoeff(sphere, excitation, n, ka)

    # derivative of associated Legendre Polynomial
    if abs(cosϑ) < 0.999999
        if n == 1
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end
    end

    if abs(sinϑ) > 1.0e-6
        t1 = An * plm[n] / sinϑ
        t2 = Bn * plm[n] / sinϑ
    end

    if cosϑ > 0.999999
        ΔSϑ = ((im)^(n - 1)) * (n * (n + 1) / 2) * (An - im * Bn) # Ruck, et. al. (3.1-12)
        ΔSϕ = ΔSϑ
    elseif cosϑ < -0.999999
        ΔSϑ = ((-im)^(n - 1)) * (n * (n + 1) / 2) * (An + im * Bn) # Ruck, et. al. (3.1-14)
        ΔSϕ = -ΔSϑ
    else
        ΔSϑ = ((im)^(n + 1)) * (t1 - im * Bn * dp) # Ruck, et. al. (3.1-6)
        ΔSϕ = ((im)^(n + 1)) * (An * dp - im * t2) # Ruck, et. al. (3.1-7)
    end

    return ΔSϑ, ΔSϕ
end
