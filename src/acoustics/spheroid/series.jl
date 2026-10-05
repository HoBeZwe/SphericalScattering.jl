
"""
    SpheroidalModes{T}

Table of the modal coefficients of a spheroidal scattering problem.

Since the spheroidal wave functions are expensive compared to the spherical ones, the coefficients are computed
once for a given scatterer and excitation and are reused for every evaluation point. Stored are the spheroidal
parameter `c`, the truncation orders `M` and `N`, the coefficients `A` of the incident expansion and the
scattering coefficients `b`, both indexed as `[m + M + 1, n + 1]` with ``-M ≤ m ≤ M`` and ``|m| ≤ n ≤ N``.
"""
struct SpheroidalModes{T}
    c::T
    M::Int
    N::Int
    A::Matrix{Complex{T}}
    b::Matrix{Complex{T}}
end


"""
    modes(sphere::Spheroid, excitation::AcousticExcitation; M::Int, N::Int, ξ=nothing, nφ=nothing, nη=nothing)

Compute the modal coefficients of the scattering problem, see [`SpheroidalModes`](@ref).

The coefficients ``A_{mn}`` of the incident expansion
``p_\\mathrm{i} = \\sum_{mn} A_{mn} R^{(1)}_{mn}(c, ξ) S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}``
are obtained by projecting the incident field, which is known in closed form, onto the angular functions on a
surface ``ξ = \\mathrm{const}``. The angular functions of equal order are orthogonal there with unit weight,
and the exponentials are orthogonal in ``φ``, so that

```math
A_{mn} = \\cfrac{1}{R^{(1)}_{mn}(c, ξ) \\, N_{mn}} \\int_{-1}^{1} \\left[ \\cfrac{1}{2π} \\int_0^{2π}
         p_\\mathrm{i} \\, \\mathrm{e}^{-\\mathrm{j}mφ} \\, \\mathrm{d}φ \\right] S_{mn}(c, η) \\, \\mathrm{d}η
```

holds with ``N_{mn} = \\int_{-1}^{1} S_{mn}^2 \\, \\mathrm{d}η``. The projection makes no assumption about the
excitation and is independent of the normalization of the angular functions, since the same normalization enters
the numerator and ``N_{mn}``.

The surface of the projection has to lie between the scatterer and the source of the incident field: beyond the
source, see [`sourceCoordinate`](@ref), the expansion in the regular wave functions does not hold. By default it
is placed at [`projectionCoordinate`](@ref), or halfway between the scatterer and the source if that is closer; a
given `ξ` reaching the source is rejected.

The truncations are determined automatically unless they are given. Since the scattering coefficients decay
once a mode is cut off at the surface, the degree is bounded by the size of the scatterer, so that
``N = \\lceil c \\, ρ \\rceil + 15`` is taken initially, where ``ρ`` is the radius of the circumscribing sphere in
units of the semifocal distance, see [`normalizedCircumradius`](@ref); ``c ρ`` is the counterpart of ``ka``. The
order `M` is not estimated but measured: the azimuthal spectrum of the incident field on the surface of the
projection is evaluated, and the orders are retained as long as they contribute more than the relative accuracy.
This matters, as every order costs a pair of calls to the spheroidal wave functions, and a nearly axial
excitation needs far fewer orders than a grazing one.

The truncation is verified rather than trusted, on the surface of the scatterer: the size of a mode of the
scattered pressure there bounds its contribution everywhere, as the outgoing radial functions decrease outward.
The modes of the last degrees have to be negligible, see [`projectedModes`](@ref). A source close to the
scatterer requires more degrees than its size suggests, as the expansion of the incident field converges the more
slowly the closer the source is. If the degree has been determined automatically, it is therefore raised by half
until the omitted modes are negligible, up to four times its initial value and as long as the projection keeps
improving: for a source very close to the scatterer the regular radial functions on the surface of the
projection become tiny, which limits the attainable accuracy. A message is printed if the relative accuracy is
not attained. A given `N`, or a `nmax` of [`Parameter`](@ref), is used as it is.

!!! note
    The degree bounds the scattered field at every radial coordinate, not only near the scatterer, since the
    scattering coefficients and the outgoing radial functions both decay with the degree. Reconstructing the
    *incident* field from the coefficients `A` is a different matter: evaluating it at a radial coordinate much
    larger than `ξ` amplifies the error of the projection, because the regular radial functions grow steeply
    with ``ξ`` at a fixed degree. For that purpose, and only for it, `ξ` should be chosen comparable to the
    largest radial coordinate of interest, and `N` large enough that ``N ≳ c ξ``.
"""
function modes(
    sphere::Spheroid,
    excitation::AcousticExcitation;
    M=nothing,
    N=nothing,
    ξ=nothing,
    nφ=nothing,
    nη=nothing,
    parameter::Parameter=Parameter(),
)

    checkExcitation(sphere, excitation)

    T = typeof(excitation.frequency)

    c = spheroidalParameter(sphere, excitation)

    eps = parameter.relativeAccuracy

    # --- the counterpart of ka bounds the degree, as the scattering coefficients decay beyond it
    adaptive = isnothing(N) && parameter.nmax < 0
    N₀ = isnothing(N) ? (parameter.nmax >= 0 ? parameter.nmax : ceil(Int, c * normalizedCircumradius(sphere)) + 15) : N

    # --- the surface of the projection, which has to lie between the scatterer and the source of the incident
    #     field: beyond the source the expansion in the regular wave functions does not hold
    ξs = sourceCoordinate(sphere, excitation)
    ξp = isnothing(ξ) ? T(min(projectionCoordinate(sphere), (sphere.ξ₀ + ξs) / 2)) : T(ξ)

    ξp < ξs || error(
        "The surface of the projection, ξ = $ξp, has to lie closer to the scatterer than the source of the incident field at ξ = $ξs.",
    )

    md, tail = projectedModes(sphere, excitation, N₀, M, ξp, nη, nφ, eps)

    # --- a source close to the scatterer requires more degrees than its size suggests: unless it is given, the
    #     degree is raised until the omitted modes are negligible, as long as the projection keeps improving.
    #     It ceases to as the regular radial functions on the surface of the projection become tiny
    saturated = false

    if adaptive
        Nlimit = 4 * N₀

        while tail > eps && md.N < Nlimit
            next, nextTail = projectedModes(sphere, excitation, min(ceil(Int, 3 * md.N / 2), Nlimit), M, ξp, nη, nφ, eps)

            if !(all(isfinite, next.A) && all(isfinite, next.b) && nextTail < tail)
                saturated = true
                break
            end

            md, tail = next, nextTail
        end
    end

    remedy = if saturated
        "more degrees do not improve the projection, the source being too close to the scatterer"
    else
        "a larger `N` may be passed explicitly"
    end

    tail > eps &&
        print("truncation may be insufficient: the modes of degree N=$(md.N) still contribute $tail on the surface; $remedy\n")

    return md
end


"""
    projectedModes(sphere::Spheroid, excitation::AcousticExcitation, Nmax::Int, M, ξp, nη, nφ, eps)

Compute the modal coefficients for the degree `Nmax` by the projection on the surface `ξp`, see [`modes`](@ref).

Returned are the [`SpheroidalModes`](@ref) and the size of the omitted modes: the L² norm of each mode of the
scattered pressure on the surface of the scatterer, ``|A_{mn} b_{mn} R^{(\\mathrm{out})}_{mn}(c, ξ_0)| \\sqrt{N_{mn}}``,
the largest among the last two degrees relative to the largest overall. Two degrees are taken, as on a disc only
modes of one parity of ``n - m`` scatter, so that the last degree alone may vanish for an axial excitation.
"""
function projectedModes(sphere::Spheroid, excitation::AcousticExcitation, Nmax::Int, M, ξp, nη, nφ, eps)

    T = typeof(excitation.frequency)

    c = spheroidalParameter(sphere, excitation)

    nη = isnothing(nη) ? 2 * Nmax + 16 : nη
    nφ = isnothing(nφ) ? 4 * Nmax + 8 : nφ               # independent of M, so that M can be measured

    ηNodes, ηWeights = gausslegendre(nη)
    φNodes = T(2π) .* (0:(nφ - 1)) ./ nφ

    R = frame(sphere)

    # --- the incident field on the surface of the projection, in the global frame
    quantity = Pressure(nothing) # the single-point methods take the quantity for dispatch only

    pᵢ = zeros(Complex{T}, nη, nφ)
    for (i, η) in enumerate(ηNodes), (j, φ) in enumerate(φNodes)
        pᵢ[i, j] = field(excitation, R * cartesianCoordinates(sphere, SVector(ξp, η, φ)), quantity)
    end

    # --- the azimuthal projections, from which the order truncation is measured
    g = [(pᵢ * cis.(-m .* φNodes)) ./ nφ for m in (-Nmax):Nmax]
    amplitude = [maximum(abs, gₘ) for gₘ in g]

    Mmax = if isnothing(M)
        significant = findall(>(eps * maximum(amplitude)), amplitude)
        maximum(abs(ind - Nmax - 1) for ind in significant)
    else
        M
    end

    Mmax <= Nmax || error("The truncation `M` of the order must not be larger than the truncation `N` of the degree.")

    A = zeros(Complex{T}, 2 * Mmax + 1, Nmax + 1)
    b = zeros(Complex{T}, 2 * Mmax + 1, Nmax + 1)

    surfaceNorm = zeros(T, 2 * Mmax + 1, Nmax + 1) # the L² norm of each mode of the scattered pressure on the surface

    M, N = Mmax, Nmax

    for m in (-M):M

        mAbs = abs(m)

        gₘ = g[m + Nmax + 1]

        # --- the angular and radial functions of this order, batched over the degree
        S = smn(mAbs, mAbs:N, c, ηNodes; spheroid=shape(sphere), normalize=false).value
        R₁ = rmn(mAbs, mAbs:N, c, [ξp]; spheroid=shape(sphere), kind=regularKind).value
        Rₒ = rmn(mAbs, mAbs:N, c, [sphere.ξ₀]; spheroid=shape(sphere), kind=outgoingKind).value

        for (k, n) in enumerate(mAbs:N)
            Sₙ = view(S, :, k)

            Nₘₙ = sum(ηWeights .* Sₙ .^ 2)

            A[m + M + 1, n + 1] = sum(ηWeights .* gₘ .* Sₙ) / (R₁[1, k] * Nₘₙ)
            b[m + M + 1, n + 1] = scatterCoeff(sphere, excitation, mAbs, n)

            surfaceNorm[m + M + 1, n + 1] = abs(A[m + M + 1, n + 1] * b[m + M + 1, n + 1] * Rₒ[1, k]) * sqrt(Nₘₙ)
        end
    end

    # --- the truncation is verified rather than trusted, on the surface of the scatterer: the outgoing radial
    #     functions decrease outward, so that the size of a mode there bounds its contribution everywhere
    tail = maximum(view(surfaceNorm, :, max(N, 1):(N + 1))) / maximum(surfaceNorm)

    return SpheroidalModes{T}(c, M, N, A, b), tail
end


"""
    sourceCoordinate(sphere::Spheroid, excitation::AcousticExcitation)

Returns the radial coordinate of the source of the incident field, beyond which its expansion in the regular
spheroidal wave functions does not hold: that of the position of a monopole, and infinity for a plane wave.
"""
sourceCoordinate(sphere::Spheroid, excitation::AcousticPlaneWave) = Inf

sourceCoordinate(sphere::Spheroid, excitation::AcousticMonopole) =
    spheroidalCoordinates(sphere, frame(sphere)' * excitation.position)[1]


"""
    seriesvalue(sphere::Spheroid, md::SpheroidalModes, point, coefficients, outgoing::Bool)

Evaluate ``\\sum_{mn} \\mathrm{coefficients}_{mn} R_{mn}(c, ξ) S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}`` at a
point given in Cartesian coordinates of the global frame.

With `outgoing=false` the regular radial functions are employed, which reproduces the incident field from the
coefficients `A`; with `outgoing=true` the outgoing ones are employed, which yields the scattered field.
"""
function seriesvalue(sphere::Spheroid, md::SpheroidalModes{T}, point, coefficients, outgoing::Bool) where {T}

    ξ, η, φ = spheroidalCoordinates(sphere, frame(sphere)' * point)

    # as for the spherical scatterers, the pressure vanishes inside. The series is evaluated regardless
    # when the regular radial functions are requested, since reproducing the incident field, which is
    # regular everywhere, is the purpose of that branch
    outgoing && ξ < sphere.ξ₀ && return Complex{T}(0.0)

    u = Complex{T}(0.0)

    for m in (-md.M):(md.M)

        mAbs = abs(m)

        S = smn(mAbs, mAbs:(md.N), md.c, [η]; spheroid=shape(sphere), normalize=false).value
        Rr = rmn(mAbs, mAbs:(md.N), md.c, [ξ]; spheroid=shape(sphere), kind=outgoing ? outgoingKind : regularKind).value

        phase = cis(m * φ)

        for (k, n) in enumerate(mAbs:(md.N))
            u += coefficients[m + md.M + 1, n + 1] * Rr[1, k] * S[1, k] * phase
        end
    end

    return u
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, quantity::Pressure; parameter::Parameter=Parameter())

Compute the pressure scattered by an oblate spheroid, for an incident acoustic excitation.

Since the surface ``ξ = ξ_0`` is a coordinate surface, the boundary condition decouples the modes, so that the
scattered pressure follows from the incident expansion by multiplying each coefficient with the scattering
coefficient of its mode and by replacing the regular radial function with the outgoing one:

```math
p_\\mathrm{s} = \\sum_{mn} A_{mn} b_{mn} R^{(\\mathrm{out})}_{mn}(c, ξ) S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}
```

The modal coefficients are passed in, see [`modes`](@ref), as the spheroidal wave functions are expensive
enough that they should be computed once per scatterer and excitation rather than per evaluation point.

The point is in Cartesian coordinates of the global frame; an arbitrary orientation of the spheroid and an
arbitrary direction or position of the excitation are accounted for.
"""
function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, quantity::Pressure; parameter::Parameter=Parameter()
)

    return seriesvalue(sphere, md, point, md.A .* md.b, true)
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Pressure; parameter::Parameter=Parameter())

Compute the pressure scattered by an oblate spheroid at all locations of `quantity`.
"""
function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Pressure; parameter::Parameter=Parameter()
)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    coefficients = md.A .* md.b

    p = progress(length(quantity.locations))

    @tasks for ind in eachindex(quantity.locations)
        F[ind] = seriesvalue(sphere, md, quantity.locations[ind], coefficients, true)
        next!(p)
    end
    finish!(p)

    return F
end


"""
    field(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Pressure; parameter::Parameter=Parameter())

Compute the total pressure in the presence of an oblate spheroid, for an incident acoustic excitation.
"""
function field(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Pressure; parameter::Parameter=Parameter()
)

    F = field(excitation, quantity; parameter=parameter)
    F .+= scatteredfield(sphere, excitation, md, quantity; parameter=parameter)

    # as for the spherical scatterers the total pressure vanishes inside. The `zeroRadius` of the
    # incident field cannot express the interior of a spheroid, hence it is masked here
    for (ind, point) in enumerate(quantity.locations)
        isinside(sphere, point) && (F[ind] = zero(eltype(F)))
    end

    return F
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::Pressure; parameter::Parameter=Parameter())

Compute the pressure scattered by an oblate spheroid, for an incident acoustic excitation.

The modal coefficients are determined automatically, see [`modes`](@ref), and are reused for all locations of
`quantity`. Pass them in explicitly in order to reuse them across several calls, as they are by far the most
expensive part of the computation.
"""
function scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::Pressure; parameter::Parameter=Parameter())

    md = modes(sphere, excitation; parameter=parameter)

    return scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
end


"""
    field(sphere::Spheroid, excitation::AcousticExcitation, quantity::Pressure; parameter::Parameter=Parameter())

Compute the total pressure in the presence of an oblate spheroid, for an incident acoustic excitation.
"""
function field(sphere::Spheroid, excitation::AcousticExcitation, quantity::Pressure; parameter::Parameter=Parameter())

    md = modes(sphere, excitation; parameter=parameter)

    return field(sphere, excitation, md, quantity; parameter=parameter)
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the far field of the pressure scattered by an oblate spheroid.

The outgoing radial functions behave as
``R^{(\\mathrm{out})}_{mn}(c, ξ) → \\mathrm{j}^{n+1} \\mathrm{e}^{-\\mathrm{j}cξ} / (cξ)`` for ``ξ → ∞``, just as
``h_n^{(2)}`` does, the limit depending neither on the order nor on ``c``. Since ``r = f \\sqrt{ξ^2 + 1 - η^2}``
tends to ``f ξ``, so that ``c ξ → k r``, the far field

```math
p^\\mathrm{sc}_\\infty = \\lim_{r → ∞} r \\, \\mathrm{e}^{\\mathrm{j} k r} p_\\mathrm{s}
    = \\cfrac{1}{k} \\sum_{mn} A_{mn} b_{mn} \\, \\mathrm{j}^{n+1} S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}
```

follows, that is, the factor ``\\mathrm{e}^{-\\mathrm{j}kr} / r`` is omitted, as for the electromagnetic
excitations. No radial function has to be evaluated.

The far field is determined by the direction of observation alone. In the limit ``η`` becomes the cosine of the
angle between the axis of the spheroid and that direction, which is why it is taken from the direction rather
than from the oblate spheroidal coordinates of the point.
"""
function scatteredfield(
    sphere::Spheroid,
    excitation::AcousticExcitation,
    md::SpheroidalModes{T},
    point,
    quantity::FarField;
    parameter::Parameter=Parameter(),
) where {T}

    coefficients = md.A .* md.b

    return farfieldvalue(sphere, md, point, coefficients) / wavenumber(excitation)
end


"""
    farfieldvalue(sphere::Spheroid, md::SpheroidalModes, point, coefficients)

Evaluate ``\\sum_{mn} \\mathrm{coefficients}_{mn} \\, \\mathrm{j}^{n+1} S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}``
for the direction of `point`, see the far field of [`scatteredfield`](@ref).
"""
function farfieldvalue(sphere::Spheroid, md::SpheroidalModes{T}, point, coefficients) where {T}

    # in the limit the angular coordinate is the cosine of the angle from the axis of the spheroid
    ~, ϑ, φ = cart2sph(frame(sphere)' * point)

    η = cos(ϑ)

    u = Complex{T}(0.0)

    for m in (-md.M):(md.M)

        mAbs = abs(m)

        S = smn(mAbs, mAbs:(md.N), md.c, [η]; spheroid=shape(sphere), normalize=false).value

        phase = cis(m * φ)

        for (k, n) in enumerate(mAbs:(md.N))
            u += coefficients[m + md.M + 1, n + 1] * im^(n + 1) * S[1, k] * phase
        end
    end

    return u
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::FarField; parameter::Parameter=Parameter())

Compute the far field of the pressure scattered by an oblate spheroid at all locations of `quantity`.
"""
function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::FarField; parameter::Parameter=Parameter()
)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    coefficients = md.A .* md.b
    k = wavenumber(excitation)

    p = progress(length(quantity.locations))

    @tasks for ind in eachindex(quantity.locations)
        F[ind] = farfieldvalue(sphere, md, quantity.locations[ind], coefficients) / k
        next!(p)
    end
    finish!(p)

    return F
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::FarField; parameter::Parameter=Parameter())

Compute the far field of the pressure scattered by an oblate spheroid, determining the modal coefficients
automatically, see [`modes`](@ref).
"""
function scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::FarField; parameter::Parameter=Parameter())

    md = modes(sphere, excitation; parameter=parameter)

    return scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
end
