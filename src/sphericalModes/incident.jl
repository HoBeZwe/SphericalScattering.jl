
"""
    field(excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field of a spherical mode.
"""
function field(excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)
    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- distinguish electric/magnetic field
    fieldType, exc = getFieldType(excitation, quantity)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(exc, point, fieldType; parameter=parameter)
    end

    return F
end



"""
    field(excitation::SphericalModeTE, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field of a TE spherical mode.
"""
function field(excitation::SphericalModeTE, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    m = excitation.m
    n = excitation.n
    k = wavenumber(excitation)

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    T = typeof(k)

    r = point_sph[1]
    ϑ = point_sph[2]
    ϕ = point_sph[3]

    kr = k * r

    Er = Complex{T}(0.0)
    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    # --- common factor
    pf = prefac(m, n)
    ZF = sqrt(ε / μ)

    # --- Hankel functions
    if excitation.c == 1        # inward
        H = hankelh1(n + T(0.5), kr) * sqrt(π / 2 / kr)
    elseif excitation.c == 2    # outward
        H = hankelh2(n + T(0.5), kr) * sqrt(π / 2 / kr)
    else
        error("Type can only be 1 or 2.")
    end

    # --- factors Eϑ and Eϕ have in common
    aux = excitation.amplitude * k * sqrt(ZF) * H * pf * cis(m * ϕ)

    # --- put things together
    Eϑ = aux * im * associatedLegendre(n, m, ϑ)
    Eϕ = -aux * derivatieAssociatedLegendre(n, m, ϑ)

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, Eϕ), point_sph) # convert to Cartesian representation
end



"""
    field(excitation::SphericalModeTE, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field of a TM spherical mode. 
"""
function field(excitation::SphericalModeTM, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    m = excitation.m
    n = excitation.n
    k = wavenumber(excitation)

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    T = typeof(k)

    r = point_sph[1]
    ϑ = point_sph[2]
    ϕ = point_sph[3]

    kr = k * r

    Er = Complex{T}(0.0)
    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    # --- common factor
    pf = prefac(m, n)
    ZF = sqrt(ε / μ)

    # --- Hankel functions
    if excitation.c == 1        # inward
        H  = hankelh1(n + T(0.5), kr)
        dH = (n + 1) * H - kr * hankelh1(n + T(1.5), kr)
    elseif excitation.c == 2    # outward
        H  = hankelh2(n + T(0.5), kr)
        dH = (n + 1) * H - kr * hankelh2(n + T(1.5), kr)
    else
        error("Type can only be 1 or 2.")
    end
    pre = sqrt(π / 2 / kr)
    H *= pre
    dH *= pre

    # --- factors Er, Eϑ, and Eϕ have in common
    mabs = abs(m)
    aux = excitation.amplitude * sqrt(ZF) * pf * cis(m * ϕ) / r

    # --- put things together
    Er = aux * n * (n + 1) * H * Plm(cos(ϑ), n, mabs; norm=Val(:normalized))
    Eϑ = aux * dH * derivatieAssociatedLegendre(n, m, ϑ)
    Eϕ = aux * dH * im * associatedLegendre(n, m, ϑ)

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, Eϕ), point_sph) # convert to Cartesian representation
end



"""
    field(excitation::SphericalModeTE, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field of a TE spherical mode.
"""
function field(excitation::SphericalModeTE, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    m = excitation.m
    n = excitation.n
    k = wavenumber(excitation)

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    T = typeof(k)

    ϑ = point_sph[2]
    ϕ = point_sph[3]


    Er = Complex{T}(0.0)
    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    # --- common factor
    pf = prefac(m, n)
    ZF = sqrt(ε / μ)

    # --- factors Eϑ and Eϕ have in common
    aux = excitation.amplitude * (im)^(n + 1) * T(1.0) * sqrt(ZF) * pf * cis(m * ϕ)

    # --- put things together
    Eϑ = aux * im * associatedLegendre(n, m, ϑ)
    Eϕ = -aux * derivatieAssociatedLegendre(n, m, ϑ)

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, Eϕ), point_sph) # convert to Cartesian representation
end



"""
    field(excitation::SphericalModeTE, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric far-field of a TM spherical mode. 
"""
function field(excitation::SphericalModeTM, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    m = excitation.m
    n = excitation.n
    k = wavenumber(excitation)

    μ = excitation.embedding.μ
    ε = excitation.embedding.ε

    T = typeof(k)

    ϑ = point_sph[2]
    ϕ = point_sph[3]

    Er = Complex{T}(0.0)
    Eϑ = Complex{T}(0.0)
    Eϕ = Complex{T}(0.0)

    # --- common factor
    pf = prefac(m, n)
    ZF = sqrt(ε / μ)

    # --- factors Er, Eϑ, and Eϕ have in common
    aux = excitation.amplitude * (im)^n * T(1.0) * sqrt(ZF) * pf * cis(m * ϕ)

    # --- put things together
    Eϑ = aux * derivatieAssociatedLegendre(n, m, ϑ)
    Eϕ = aux * im * associatedLegendre(n, m, ϑ)

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, Eϕ), point_sph) # convert to Cartesian representation
end



"""
    prefac(m::T, n::T)

Prefactor for both types of spherical vector wave functions.

Note the factor (-m / abs(m))^m is now part of the associated Legendre polynomials via the Condon-Shortley phase.
"""
function prefac(m::T, n::T) where {T<:Integer}

    return 1.0 / sqrt(n * (n + 1) * 2π)
end



"""
    associatedLegendre(n::T, m::T, ϑ::F)

Compute the normalized associated Legendre polynomials according to the definition by Hansen (A1.25) times m divided by sin(ϑ), 
but including the Condon-Shortley phase.
By using the recursion relations (A1.34a) all corner cases are accounted for correctly.

m * bar(P)_n_|m|(cos(ϑ)) / sin(ϑ)

ϑ ∈ [0, π] assumed.
"""
function associatedLegendre(n::T, m::T, ϑ::F) where {T<:Integer,F<:Real}

    m == 0 && return F(0.0)

    cosϑ = cos(ϑ)
    sinϑ = sin(ϑ)
    mabs = abs(m)
    signm = sign(m)

    mabs == n && return signm * -F(0.5) * cosϑ * (sqrt((n - mabs + 1) * (n + mabs)) * Plm(cosϑ, n, mabs - 1; norm=Val(:normalized))) +
           signm * mabs * sinϑ * Plm(cosϑ, n, mabs; norm=Val(:normalized))

    return signm *
           -F(0.5) *
           cosϑ *
           (
               sqrt((n - mabs + 1) * (n + mabs)) * Plm(cosϑ, n, mabs - 1; norm=Val(:normalized)) +
               sqrt((n - mabs) * (n + mabs + 1)) * Plm(cosϑ, n, mabs + 1; norm=Val(:normalized))
           ) + signm * mabs * sinϑ * Plm(cosϑ, n, mabs; norm=Val(:normalized))
end



"""
    derivatieAssociatedLegendre(n::T, m::T, ϑ::F) where {T<:Integer,F<:Real}

Compute the derivatives of the normalized associated Legendre polynomials according to the definition by Hansen (A1.25), 
but including the Condon-Shortley phase.
By using the recursion relations (A1.34b) all corner cases are accounted for correctly.

d bar(P)_n_|m|(cos(ϑ)) / dϑ

ϑ ∈ [0, π] assumed.
"""
function derivatieAssociatedLegendre(n::T, m::T, ϑ::F) where {T<:Integer,F<:Real}

    cosϑ = cos(ϑ)
    mabs = abs(m)

    mabs == 0 && return sqrt((n) * (n + 1)) * Plm(cosϑ, n, 1; norm=Val(:normalized))

    mabs == n && return -F(0.5) * (sqrt((n - mabs + 1) * (n + mabs)) * Plm(cosϑ, n, mabs - 1; norm=Val(:normalized)))

    return -F(0.5) * (
        sqrt((n - mabs + 1) * (n + mabs)) * Plm(cosϑ, n, mabs - 1; norm=Val(:normalized)) -
        sqrt((n - mabs) * (n + mabs + 1)) * Plm(cosϑ, n, mabs + 1; norm=Val(:normalized))
    )
end
