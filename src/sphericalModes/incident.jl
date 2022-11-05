
"""
    field(excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field of a spherical mode.
"""
function field(excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

    T = typeof(excitation.wavenumber)
    F = zeros(SVector{3,Complex{T}}, length(quantity.locations))

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
    k = excitation.wavenumber

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

    # --- Hankel functions
    if excitation.c == 1        # inward
        H = hankelh1(n + T(0.5), kr) * sqrt(π / 2 / kr)
    elseif excitation.c == 2    # outward
        H = hankelh2(n + T(0.5), kr) * sqrt(π / 2 / kr)
    else
        error("Type can only be 1 or 2.")
    end

    # --- factors Eϑ and Eϕ have in common
    mabs = abs(m)
    aux =
        excitation.amplitude *
        k *
        sqrt(sqrt(μ / ε)) *
        H *
        pf *
        exp(-im * m * ϕ) *
        sqrt((2 * n + 1) / 2 * factorial(n - mabs) / factorial(n + mabs))

    # --- put things together
    Eϑ = aux * im * associatedLegendre(n, m, ϑ)
    Eϕ = aux * derivatieAssociatedLegendre(n, m, ϑ)

    #return SVector(Er, Eϑ, Eϕ)
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
    k = excitation.wavenumber

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
    aux =
        excitation.amplitude *
        sqrt(sqrt(μ / ε)) *
        pf *
        exp(-im * m * ϕ) *
        sqrt((2 * n + 1) / 2 * factorial(n - mabs) / factorial(n + mabs)) / r

    # --- put things together
    Er = aux * n * (n + 1) * H * (-1)^mabs * Plm(cos(ϑ), n, mabs)
    Eϑ = aux * dH * derivatieAssociatedLegendre(n, m, ϑ)
    Eϕ = -aux * dH * im * associatedLegendre(n, m, ϑ)

    #return SVector(Er, Eϑ, Eϕ)
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
    k = excitation.wavenumber

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

    # --- factors Eϑ and Eϕ have in common
    mabs = abs(m)
    aux =
        -excitation.amplitude *
        (im)^(n + 1) *
        T(1.0) *
        sqrt(sqrt(μ / ε)) *
        pf *
        exp(-im * m * ϕ) *
        sqrt((2 * n + 1) / 2 * factorial(n - mabs) / factorial(n + mabs))

    # --- put things together
    Eϑ = aux * im * associatedLegendre(n, m, ϑ)
    Eϕ = aux * derivatieAssociatedLegendre(n, m, ϑ)

    #return SVector(Er, Eϑ, Eϕ)
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
    k = excitation.wavenumber

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

    # --- factors Er, Eϑ, and Eϕ have in common
    mabs = abs(m)
    aux =
        excitation.amplitude *
        (-im)^n *
        T(1.0) *
        sqrt(sqrt(μ / ε)) *
        pf *
        exp(-im * m * ϕ) *
        sqrt((2 * n + 1) / 2 * factorial(n - mabs) / factorial(n + mabs))

    # --- put things together
    Eϑ = -aux * derivatieAssociatedLegendre(n, m, ϑ)
    Eϕ = aux * im * associatedLegendre(n, m, ϑ)

    #return SVector(Er, Eϑ, Eϕ)
    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Er, Eϑ, Eϕ), point_sph) # convert to Cartesian representation
end



"""
    prefac(m::T, n::T)

Prefactor for both types of spherical vector wave functions
"""
function prefac(m::T, n::T) where {T<:Integer}

    aux = 1.0 / sqrt(n * (n + 1) * 2π)

    # special case m = 0
    m == 0 && return aux

    # all other cases
    return aux * (-m / abs(m))^m
end



"""
    associatedLegendre(n::T, m::T, ϑ::F)

Compute the associated Legendre polynomials according to the definition by Hansen (A1.25) times m divided by sin(ϑ):

m * P_n_|m|(cos(ϑ)) / sin(ϑ)

The special values for ϑ ∈ {0, π/2, π} are treated properly.

ϑ ∈ [0, π] assumed
"""
function associatedLegendre(n::T, m::T, ϑ::F) where {T<:Integer,F<:Real}

    cosϑ = cos(ϑ)
    sinϑ = sin(ϑ)
    mabs = abs(m)

    # --- handle limit cases
    if isapprox(ϑ, 0.0; atol=1e-6)

        mabs == 1 && return (-1)^mabs * F(m * n * (n + 1) / 2)
        return F(0.0)
    end

    if isapprox(ϑ, π; rtol=1e-6)

        mabs == 1 && return (-1)^mabs * F(m * (-1)^(n + 1) * n * (n + 1) / 2)
        return F(0.0)
    end

    return (-1)^mabs * m * Plm(cosϑ, n, mabs) / sinϑ
end

# function associatedLegendre(n::T, m::T, ϑ::F) where {T <: Integer, F <: Real}

#     cosϑ = cos(ϑ)
#     sinϑ = sin(ϑ)
#     mabs = abs(m)

#     mabs == 0 && return F(0.0)

#     return (-1)^mabs * m * sinϑ * Plm(cosϑ, n, mabs) + 0.5 * cosϑ * ((n - mabs + 1)*(n + mabs) * (-1)^(mabs-1) * Plm(cosϑ, n, mabs-1) + (-1)^(mabs+1) * Plm(cosϑ, n, mabs+1))
# end


"""

d P_n_|m|(cos(ϑ)) / dϑ
"""
function derivatieAssociatedLegendre(n::T, m::T, ϑ::F) where {T<:Integer,F<:Real}

    cosϑ = cos(ϑ)
    mabs = abs(m)

    mabs == 0 && return Plm(cosϑ, n, 1)
    m == n && return F(0.5) * ((n - mabs + 1) * (n + mabs) * (-1)^(mabs - 1) * Plm(cosϑ, n, mabs - 1))
    return F(0.5) * ((n - mabs + 1) * (n + mabs) * (-1)^(mabs - 1) * Plm(cosϑ, n, mabs - 1) - (-1)^(mabs + 1) * Plm(cosϑ, n, mabs + 1))
end
