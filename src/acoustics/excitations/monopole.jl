
struct AcousticMonopole{T,R,C} <: AcousticExcitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
end


"""
    field(excitation::AcousticMonopole, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure of an acoustic monopole, or its Dirichlet trace.

The Dirichlet trace is the pressure itself. Its locations are taken as they are given and are assumed to lie on
the surface of interest, whereas the traces of the scattered and of the total field are evaluated on the surface
of the sphere.
"""
function field(excitation::AcousticMonopole, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter(), zeroRadius=0.0)

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
    field(excitation::AcousticMonopole, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace of the pressure of an acoustic monopole for the normals of `quantity`.

The locations are taken as they are given and are assumed to lie on the surface of interest, whereas the traces
of the scattered and of the total field are evaluated on the surface of the sphere.
"""
function field(excitation::AcousticMonopole, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute trace
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity.normals[ind], quantity; parameter=parameter)
    end

    return F
end


"""
    field(excitation::AcousticMonopole, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure ``p_\\mathrm{i} = A \\mathrm{e}^{-\\mathrm{j}kR} / (4πR)`` of an acoustic monopole, which is
at the same time its Dirichlet trace ``γ_0 p_\\mathrm{i}``.

Here, ``R`` denotes the distance to the monopole, so that the pressure is the free-space Green's function of the
Helmholtz equation scaled by the amplitude ``A``. The field is singular at the position of the monopole.

The point is in Cartesian coordinates.
"""
function field(excitation::AcousticMonopole, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    R = norm(point - excitation.position)

    return a / (4 * π) * cis(-k * R) / R
end


"""
    field(excitation::AcousticMonopole, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace ``γ_1 p_\\mathrm{i} = \\hat{n} ⋅ ∇ p_\\mathrm{i}`` of the pressure of an acoustic
monopole for the given `normal`.

Since the pressure depends on the position solely via the distance ``R`` to the monopole,

```math
∇ p_\\mathrm{i} = \\frac{∂ p_\\mathrm{i}}{∂R} \\hat{R}
                = -\\frac{A}{4π} (1 + \\mathrm{j}kR) \\frac{\\mathrm{e}^{-\\mathrm{j}kR}}{R^2} \\hat{R}
```

holds, where ``\\hat{R}`` points from the monopole to the point. For ``kR ≫ 1`` this approaches
``-\\mathrm{j} k (\\hat{R} ⋅ \\hat{n}) p_\\mathrm{i}``, the local plane-wave result.

The point and the normal are in Cartesian coordinates, the latter being a unit vector.
"""
function field(excitation::AcousticMonopole, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    d = point - excitation.position # vector pointing from the monopole to the point
    R = norm(d)

    return -a / (4 * π) * (1 + im * k * R) * cis(-k * R) / R^2 * dot(normal, d / R)
end


"""
    field(excitation::AcousticMonopole, quantity::FarField; parameter::Parameter=Parameter())

Compute the far field of an acoustic monopole.

In contrast to a plane wave, a monopole does possess a far field. It is determined by the direction of
observation alone, hence no locations are suppressed and a `zeroRadius` is not taken into account.
"""
function field(excitation::AcousticMonopole, quantity::FarField; parameter::Parameter=Parameter(), zeroRadius=0.0)

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute field
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity; parameter=parameter)
    end

    return F
end


"""
    field(excitation::AcousticMonopole, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the far field ``A \\mathrm{e}^{\\mathrm{j} k \\hat{r} ⋅ \\mathbf{r}_0} / (4π)`` of an acoustic monopole
located at ``\\mathbf{r}_0``.

Since ``|\\mathbf{r} - \\mathbf{r}_0| = r - \\hat{r} ⋅ \\mathbf{r}_0 + 𝒪(1/r)``, the far field follows from the
pressure by dropping the factor ``\\mathrm{e}^{-\\mathrm{j}kr} / r``, that is, the returned quantity is
``\\lim_{r → ∞} r \\, \\mathrm{e}^{\\mathrm{j}kr} p_\\mathrm{i}``, as for the electromagnetic excitations. Its
magnitude is the same in all directions, as it has to be for a point source: the position of the monopole enters
the phase alone.

The point is in Cartesian coordinates, only its direction being relevant.
"""
function field(excitation::AcousticMonopole, point, quantity::FarField; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    r̂ = normalize(point) # only the direction of observation is relevant

    return a / (4 * π) * cis(k * dot(r̂, excitation.position))
end


"""
    symmetryAxis(excitation::AcousticMonopole)

Returns the direction towards the monopole, about which the scattered field is rotationally symmetric.
"""
symmetryAxis(excitation::AcousticMonopole) = normalize(excitation.position)


"""
    incidentCoeff(excitation::AcousticMonopole, n::Int)

Compute the coefficient ``e_n = -\\mathrm{j} k \\, h_n^{(2)}(k R_0) / (4π)`` of the n-th term of the incident
expansion, where ``R_0`` denotes the distance of the monopole from the center of the sphere.

It follows from the addition theorem of the free-space Green's function,

```math
\\frac{\\mathrm{e}^{-\\mathrm{j}k|\\mathbf{r} - \\mathbf{r}_0|}}{4π|\\mathbf{r} - \\mathbf{r}_0|}
    = \\frac{-\\mathrm{j}k}{4π} \\sum_n (2n+1) j_n(k r_<) h_n^{(2)}(k r_>) P_n(\\cos\\vartheta)
```

with ``r_< = \\min(r, R_0)``, ``r_> = \\max(r, R_0)`` and ``\\vartheta`` measured from the monopole. Since the
monopole lies outside the sphere, ``r_> = R_0`` holds at its surface, which is where the boundary condition
determines the scattering coefficients.
"""
function incidentCoeff(excitation::AcousticMonopole, n::Int)

    T = typeof(excitation.frequency)

    k = wavenumber(excitation)
    kR₀ = k * norm(excitation.position)

    s = sqrt(π / 2 / kR₀)

    return -im * k / (4 * π) * s * hankelh2(n + T(0.5), kR₀) # spherical Hankel function
end


"""
    checkExcitation(scatterer::Scatterer{<:AcousticBoundary}, excitation::AcousticMonopole)

Ensure that the monopole lies outside the scatterer, as the expansion of its field assumes.
"""
function checkExcitation(scatterer::Scatterer{<:AcousticBoundary}, excitation::AcousticMonopole)

    isinside(scatterer, excitation.position) &&
        error("The monopole has to be located outside the scatterer, as the expansion of its field assumes.")

    return nothing
end



"""
    incidentCoefficients(sphere::Spheroid, excitation::AcousticMonopole, N::Int)

Compute the coefficients of the expansion of the field of the monopole in the spheroidal wave functions, see
[`incidentCoefficients`](@ref). With the position ``(ξ_s, η_s, φ_s)`` of the monopole in the frame of the
spheroid, the addition theorem of the free-space Green's function (Flammer, 1957) yields

```math
A_{mn} = -\\cfrac{\\mathrm{j} k a}{2π} \\, \\cfrac{S_{|m|n}(c, η_s) \\, R^{(\\mathrm{out})}_{|m|n}(c, ξ_s)}{N_{|m|n}}
         \\, \\mathrm{e}^{-\\mathrm{j}mφ_s} \\,,
```

valid for ``ξ < ξ_s``, which includes the surface of the scatterer, as the monopole lies outside. It is the
counterpart of ``-\\mathrm{j}k / (4π) \\, (2n+1) \\, h_n^{(2)}(k r_0)`` for a sphere.
"""
function incidentCoefficients(sphere::Spheroid, excitation::AcousticMonopole, N::Int)

    T = typeof(excitation.frequency)

    c = spheroidalParameter(sphere, excitation)
    k = wavenumber(excitation)

    # --- the position of the monopole in the spheroidal coordinates
    ξ, η, φ = spheroidalCoordinates(sphere, frame(sphere)' * excitation.position)

    A = zeros(Complex{T}, 2 * N + 1, N + 1)

    for mAbs in 0:N
        S = smn(mAbs, mAbs:N, c, [η]; spheroid=shape(sphere), normalize=false).value

        all(iszero, S) && continue # a monopole on the axis excites the order zero alone

        Rₒ = rmn(mAbs, mAbs:N, c, [ξ]; spheroid=shape(sphere), kind=outgoingKind).value

        for (kₙ, n) in enumerate(mAbs:N), m in (iszero(mAbs) ? (0,) : (mAbs, -mAbs))
            A[m + N + 1, n + 1] = -im * k * excitation.amplitude / (2π) * S[1, kₙ] * Rₒ[1, kₙ] / angularNorm(mAbs, n) * cis(-m * φ)
        end
    end

    return A
end
