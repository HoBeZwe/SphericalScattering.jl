
"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())
    
Compute the electric field scattered by a PEC sphere, for an incident plane wave.
"""
function scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)
    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- rotate coordinates
    points = rotate(excitation, quantity.locations; inverse=true)

    p = progress(length(points))

    # --- compute field in Cartesian representation
    @tasks for ind in eachindex(points)
        F[ind] = scatteredfield(sphere, excitation, points[ind], quantity; parameter=parameter)
        next!(p)
    end
    finish!(p)

    # --- rotate resulting field
    rotate!(excitation, F; inverse=false)

    return F
end



"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, point, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC or dielectric sphere, for an incident plane wave
travelling in +z-direction with E-field polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::Sphere, excitation::PlaneWave, point, quantity::Field; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]
    r = point_sph[1]

    T = typeof(excitation.frequency)

    eps = parameter.relativeAccuracy

    ST = SVector{3,Complex{T}}
    F = ST(0.0, 0.0, 0.0)

    A₀ = amplitude(sphere, excitation::PlaneWave, quantity, r)

    A₀ == 0.0 && return F # Inside of PEC return zero field

    δF = T(Inf)
    n = 0

    sinϑ = abs(sin(point_sph[2]))  # note: theta only defined from from 0 to pi
    cosϑ = cos(point_sph[2])       # ok for theta > pi
    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    # first two values of the Associated Legendre Polynomial
    plm = Vector{T}()
    push!(plm, -sinϑ)
    push!(plm, -T(3.0) * sinϑ * cosϑ)

    try
        while δF > eps || n < 10
            n += 1

            coeffs = scatterCoeff(sphere, excitation, n)
            expansions = expansion(sphere, excitation, quantity, r, plm, cosϑ, sinϑ, n)
            ΔF = Δfieldₙ(sphere, excitation, quantity, r, coeffs, expansions, sinϕ, cosϕ, n)

            F += ΔF

            δF = norm(ΔF) / norm(F) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    return convertSpherical2Cartesian(A₀ .* F, point_sph)
end



function Δfieldₙ(sphere, excitation::PlaneWave, quantity::ElectricField, r, coeffs, expansion, sinϕ, cosϕ, n)

    (Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ) = expansion

    k = wavenumber(sphere, excitation, r)
    kr = k * r

    aₙ, bₙ = scatterCoeff_of_layer(sphere, r, coeffs)

    ΔEr = +(cosϕ / (im * kr^2)) * aₙ * Nn_r
    ΔEϑ = -(cosϕ / kr) * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
    ΔEϕ = +(sinϕ / kr) * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)

    return SVector(ΔEr, ΔEϑ, ΔEϕ)
end



function Δfieldₙ(sphere, excitation::PlaneWave, quantity::MagneticField, r, coeffs, expansion, sinϕ, cosϕ, n)

    (Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ) = expansion

    k = wavenumber(sphere, excitation, r)
    kr = k * r

    aₙ, bₙ = scatterCoeff_of_layer(sphere, r, coeffs)

    ΔHr = +(sinϕ / (im * kr^2)) * bₙ * Nn_r
    ΔHϑ = -(sinϕ / kr) * (aₙ * Mn_ϑ + bₙ * Nn_ϑ)
    ΔHϕ = -(cosϕ / kr) * (aₙ * Mn_ϕ + bₙ * Nn_ϕ)

    return SVector(ΔHr, ΔHϑ, ΔHϕ)
end



function Δfieldₙ(sphere, excitation::PlaneWave, quantity::FarField, r, coeffs, expansion, sinϕ, cosϕ, n)

    (~, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ) = expansion

    k = wavenumber(excitation)
    aₙ = coeffs[1]
    bₙ = coeffs[2]

    # See Jin, (7.4.44)-(7.4.45)
    # We deviate from Jin by replacing exp(-im*kr) / kr with 1/k
    ΔEϑ = -im * cosϕ * 1 / k * im^n * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
    ΔEϕ = +im * sinϕ * 1 / k * im^n * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)

    return SVector(eltype(ΔEϑ)(0.0), ΔEϑ, ΔEϕ)
end



"""
    amplitude(sphere, excitation::PlaneWave, quantity::MagneticField, r)

Returns ``H₀/ηᵢ``, where ``H₀`` is the magnetic field of the incident plane wave.
For PEC layers, it returns ``0``.
"""
function amplitude(sphere, excitation::PlaneWave, quantity::MagneticField, r)
    η = impedance(sphere, excitation, r)

    if η == 0.0 # PEC case
        return η * excitation.amplitude # return zero of correct type
    else
        return 1 / η * excitation.amplitude
    end
end



"""
    amplitude(sphere, excitation::PlaneWave, quantity::MagneticField, r)

Returns ``E₀``, where ``E₀`` is the electric field of the incident plane wave.
For PEC layers, it returns ``0``.
"""
function amplitude(sphere, excitation::PlaneWave, quantity::ElectricField, r)
    η = impedance(sphere, excitation, r)

    if η == 0.0 # PEC case
        return η * excitation.amplitude # return zero of correct type
    else
        return excitation.amplitude
    end
end


"""
    amplitude(sphere, excitation::PlaneWave, quantity::MagneticField, r)

Returns ``E₀``, where ``E₀`` is the electric field of the incident plane wave.
"""
function amplitude(sphere, excitation::PlaneWave, quantity::FarField, r)
    return excitation.amplitude
end



"""
    scatterCoeff_of_layer(sphere, r, coeffs)

Given a tuple of scattering coefficients `coeffs`, it returns the two scattering
coefficients of the layer to which `r` belongs.
"""
function scatterCoeff_of_layer(sphere, r, coeffs)
    cmp_index = layer(sphere, r)
    aₙ = coeffs[end - 2 * (cmp_index - 1) - 1]
    bₙ = coeffs[end - 2 * (cmp_index - 1)]
    return aₙ, bₙ
end



"""
    scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int)

Compute scattering coefficients for a plane wave travelling in +z-direction 
with polarization in x-direction.
"""
function scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int)

    T = typeof(excitation.frequency)

    k = wavenumber(excitation)

    ka = k * sphere.radius

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

    # PEC: ensure field is zero inside
    # We could enforce this somewhat differently
    # But by ensuring that the number of coefficients = 2 * number of layers
    # the code is more elegant

    cₙ = T(0) * aₙ
    dₙ = T(0) * aₙ

    return aₙ, bₙ, cₙ, dₙ
end



function scatterCoeff(sphere::DielectricSphere, excitation::PlaneWave, n::Int)

    f = excitation.frequency
    T = typeof(f)

    ε2 = excitation.embedding.ε
    μ2 = excitation.embedding.μ

    ε1 = sphere.filling.ε
    μ1 = sphere.filling.μ

    c2 = 1 / sqrt(ε2 * μ2)
    c1 = 1 / sqrt(ε1 * μ1)

    k2 = 2π * f / c2
    k1 = 2π * f / c1

    εᵣ = ε1 / ε2
    μᵣ = μ1 / μ2

    k₂a = k2 * sphere.radius
    s₂ = sqrt(π / 2 / k₂a)

    Ĵ₂  = k₂a * s₂ * besselj(n + T(0.5), k₂a)   # Riccati-Bessel function
    Ĥ₂  = k₂a * s₂ * hankelh2(n + T(0.5), k₂a)  # Riccati-Hankel function
    Ĵ₂2 = k₂a * s₂ * besselj(n - T(0.5), k₂a)   # for derivate needed
    Ĥ₂2 = k₂a * s₂ * hankelh2(n - T(0.5), k₂a)  # for derivate needed

    # Use recurrence relationship
    dĴ₂ = (Ĵ₂2 - n / k₂a * Ĵ₂)    # derivative Riccati-Bessel function
    dĤ₂ = (Ĥ₂2 - n / k₂a * Ĥ₂)    # derivative Riccati-Hankel function

    k₁a = k1 * sphere.radius
    s₁ = sqrt(π / 2 / k₁a)

    Ĵ₁ = k₁a * s₁ * besselj(n + T(0.5), k₁a)   # Riccati-Bessel function
    #Ĥ₁  = k₁a * s₁ * hankelh2(n + T(0.5), k₁a)  # Riccati-Hankel function
    Ĵ₁2 = k₁a * s₁ * besselj(n - T(0.5), k₁a)   # for derivate needed
    #Ĥ₁2 = k₁a * s₁ * hankelh2(n - T(0.5), k₁a)  # for derivate needed

    # Use recurrence relationship
    dĴ₁ = (Ĵ₁2 - n / k₁a * Ĵ₁)    # derivative Riccati-Bessel function
    #dĤ₁ = (Ĥ₁2 - n / k₁a * Ĥ₁)    # derivative Riccati-Hankel function

    pF = im^(-T(n)) * (2 * n + 1) / (n * (n + 1))

    aₙ = pF * (√εᵣ * dĴ₂ * Ĵ₁ - √μᵣ * Ĵ₂ * dĴ₁) / (√μᵣ * Ĥ₂ * dĴ₁ - √εᵣ * dĤ₂ * Ĵ₁)     # Jin (7.4.65)
    bₙ = pF * (√μᵣ * dĴ₂ * Ĵ₁ - √εᵣ * Ĵ₂ * dĴ₁) / (√εᵣ * Ĥ₂ * dĴ₁ - √μᵣ * dĤ₂ * Ĵ₁)     # Jin (7.4.66)
    cₙ = pF * (im * √εᵣ * μᵣ) / (√μᵣ * Ĥ₂ * dĴ₁ - √εᵣ * dĤ₂ * Ĵ₁)      # Jin (7.4.66)
    dₙ = pF * (im * √εᵣ * μᵣ) / (√εᵣ * Ĥ₂ * dĴ₁ - √μᵣ * dĤ₂ * Ĵ₁)      # Jin (7.4.66)

    return aₙ, bₙ, cₙ, dₙ
end



"""
    expansion(sphere::Sphere, excitation::PlaneWave, quantity::Field, r, plm, cosϑ, sinϑ, n::Int) 

Compute functional dependencies of the Mie series for a plane wave
travelling in +z-direction with polarization in x-direction.
"""
function expansion(sphere::Sphere, excitation::PlaneWave, quantity::Field, r, plm, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.frequency)

    if quantity isa FarField
        B = T(1.0)
        dB = T(1.0)
    else
        k = wavenumber(sphere, excitation, r)
        kr = k * r
        s = sqrt(π / 2 / kr)

        if r >= sphere.radius
            B  = kr * s * hankelh2(n + T(0.5), kr)     # Riccati-Hankel function
            B2 = kr * s * hankelh2(n - T(0.5), kr)
        else
            B  = kr * s * besselj(n + T(0.5), kr)     # Riccati-Hankel function
            B2 = kr * s * besselj(n - T(0.5), kr)
        end

        dB = B2 - n / kr * B          # derivative of Riccati-Hankel function
    end

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = B * p / sinϑ
        Mn_ϕ = B * dp

        Nn_ϑ = dB * dp
        Nn_ϕ = dB * p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -B * aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -dB * aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = B * aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -dB * aux
        Nn_ϕ = -Nn_ϑ
    end

    if !(quantity isa FarField)
        Nn_r = n * (n + 1) * B * p
        Nn_ϑ = im * Nn_ϑ
        Nn_ϕ = im * Nn_ϕ
    else
        Nn_r = T(0)
    end

    return Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end
