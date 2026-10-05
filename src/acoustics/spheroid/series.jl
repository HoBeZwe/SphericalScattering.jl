
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
    modes(sphere::Spheroid, excitation::AcousticExcitation; M=nothing, N=nothing, parameter::Parameter=Parameter())

Compute the modal coefficients of the scattering problem, see [`SpheroidalModes`](@ref).

The incident field is expanded as
``p_\\mathrm{i} = \\sum_{mn} A_{mn} R^{(1)}_{mn}(c, ξ) S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}``, the coefficients
``A_{mn}`` being known analytically for the plane wave and the monopole, see [`incidentCoefficients`](@ref). Since
the surface ``ξ = ξ_0`` is a coordinate surface, the boundary condition decouples the modes, each being scattered
with the coefficient ``b_{mn}`` of [`scatterCoeff`](@ref).

The truncations are determined automatically unless they are given. Since the scattering coefficients decay
once a mode is cut off at the surface, the degree is bounded by the size of the scatterer, so that
``N = \\lceil c \\, ρ \\rceil + 15`` is taken initially, where ``ρ`` is the radius of the circumscribing sphere in
units of the semifocal distance, see [`normalizedCircumradius`](@ref); ``c ρ`` is the counterpart of ``ka``. The
order `M` is not estimated but measured, see [`modesForDegree`](@ref): a nearly axial excitation needs far fewer
orders than a grazing one, which matters, as every order costs calls to the spheroidal wave functions.

The truncation is verified rather than trusted, on the surface of the scatterer: the size of a mode of the
scattered pressure there bounds its contribution everywhere, as the outgoing radial functions decrease outward.
The modes of the last degrees have to be negligible. A source close to the scatterer requires more degrees than
its size suggests, as the expansion of its field converges the more slowly the closer the source is. If the
degree has been determined automatically, it is therefore raised by half until the omitted modes are negligible,
up to four times its initial value and as long as the result keeps improving: for a source very close to the
scatterer the radial functions of high degree eventually exceed the range of the floating-point numbers, which
limits the attainable accuracy. A message is printed if the relative accuracy is not attained. A given `N`, or a
`nmax` of [`Parameter`](@ref), is used as it is.
"""
function modes(sphere::Spheroid, excitation::AcousticExcitation; M=nothing, N=nothing, parameter::Parameter=Parameter())

    checkExcitation(sphere, excitation)

    c = spheroidalParameter(sphere, excitation)

    eps = parameter.relativeAccuracy

    # --- the counterpart of ka bounds the degree, as the scattering coefficients decay beyond it
    adaptive = isnothing(N) && parameter.nmax < 0
    N₀ = isnothing(N) ? (parameter.nmax >= 0 ? parameter.nmax : ceil(Int, c * normalizedCircumradius(sphere)) + 15) : N

    md, tail = modesForDegree(sphere, excitation, N₀, M, eps)

    isnothing(md) && error(
        "The modal coefficients up to the degree N=$N₀ exceed the range of the floating-point numbers, as the wave functions of high degree and order do, in particular for a source close to the scatterer: choose a smaller degree.",
    )

    # --- a source close to the scatterer requires more degrees than its size suggests: unless it is given, the
    #     degree is raised until the omitted modes are negligible, as long as the result keeps improving. It
    #     ceases to as the wave functions of high degree exceed the range of the floating-point numbers
    saturated = false

    if adaptive
        Nlimit = 4 * N₀

        while tail > eps && md.N < Nlimit
            next, nextTail = modesForDegree(sphere, excitation, min(ceil(Int, 3 * md.N / 2), Nlimit), M, eps)

            if isnothing(next) || !(nextTail < tail)
                saturated = true
                break
            end

            md, tail = next, nextTail
        end
    end

    remedy = if saturated
        "more degrees do not improve the result, the wave functions of high degree exceeding the range of the floating-point numbers"
    else
        "a larger `N` may be passed explicitly"
    end

    tail > eps &&
        print("truncation may be insufficient: the modes of degree N=$(md.N) still contribute $tail on the surface; $remedy\n")

    return md
end


"""
    modesForDegree(sphere::Spheroid, excitation::AcousticExcitation, Nmax::Int, M, eps)

Compute the modal coefficients for the degree `Nmax`, see [`modes`](@ref).

The incident coefficients and the scattering coefficients are computed order by order, the radial functions
batched over the degree and shared by the orders ``±m``; orders the excitation does not excite, such as all but the
order zero for an axial excitation, are skipped. Unless `M` is given, the orders are computed in increasing ``|m|``
until two consecutive ones are negligible, measured by the size of their modes on the surface of the scatterer as
below, and truncated at the last one that is not. The azimuthal content of the incident field decays beyond its
angular bandwidth, so that the orders of high ``|m|`` are not computed at all: this saves their cost, and it keeps
clear of the wave functions of high order, which exceed the range of the floating-point numbers once ``n + m``
approaches about 170.

Returned are the [`SpheroidalModes`](@ref) and the size of the omitted modes: the L² norm of each mode of the
scattered pressure on the surface of the scatterer, ``|A_{mn} b_{mn} R^{(\\mathrm{out})}_{mn}(c, ξ_0)| \\sqrt{N_{mn}}``,
the largest among the last two degrees relative to the largest overall. Two degrees are taken, as on a disc only
modes of one parity of ``n - m`` scatter, so that the last degree alone may vanish for an axial excitation.

If any of the retained coefficients is not finite, as the wave functions of high degree exceed the range of the
floating-point numbers for a source very close to the scatterer, `(nothing, Inf)` is returned instead.
"""
function modesForDegree(sphere::Spheroid, excitation::AcousticExcitation, Nmax::Int, M, eps)

    T = typeof(excitation.frequency)

    c = spheroidalParameter(sphere, excitation)

    isnothing(M) || M <= Nmax || error("The truncation `M` of the order must not be larger than the truncation `N` of the degree.")

    Aall = incidentCoefficients(sphere, excitation, Nmax)

    ball = zeros(Complex{T}, 2 * Nmax + 1, Nmax + 1)
    sizeAll = zeros(T, 2 * Nmax + 1, Nmax + 1) # the L² norm of each mode of the scattered pressure on the surface

    orderSize = zeros(T, Nmax + 1) # the size of the largest mode of each order, indexed by |m| + 1
    negligible = 0

    for mAbs in 0:(isnothing(M) ? Nmax : M)

        rows = iszero(mAbs) ? (Nmax + 1,) : (Nmax + 1 + mAbs, Nmax + 1 - mAbs)

        if !all(row -> all(iszero, view(Aall, row, :)), rows) # skip an order the excitation does not excite

            # --- the radial functions on the surface, batched over the degree and shared by the orders ±m
            R₁ = rmn(mAbs, mAbs:Nmax, c, [sphere.ξ₀]; spheroid=shape(sphere), kind=regularKind)
            Rₒ = rmn(mAbs, mAbs:Nmax, c, [sphere.ξ₀]; spheroid=shape(sphere), kind=outgoingKind)

            for (k, n) in enumerate(mAbs:Nmax)
                bₙ = boundaryRatio(sphere, R₁.value[1, k], R₁.derivative[1, k], Rₒ.value[1, k], Rₒ.derivative[1, k])

                for row in rows
                    ball[row, n + 1] = bₙ
                    sizeAll[row, n + 1] = abs(Aall[row, n + 1] * bₙ * Rₒ.value[1, k]) * sqrt(angularNorm(mAbs, n))
                end
            end

            orderSize[mAbs + 1] = maximum(row -> maximum(view(sizeAll, row, :)), rows)

            # a coefficient beyond the range of the floating-point numbers, or one formed from such values
            all(row -> all(isfinite, view(Aall, row, :)) && all(isfinite, view(ball, row, :)), rows) &&
            isfinite(orderSize[mAbs + 1]) || return nothing, T(Inf)
        end

        # --- the orders are measured: once two consecutive ones are negligible, the azimuthal content has decayed
        isnothing(M) || continue
        negligible = orderSize[mAbs + 1] > eps * maximum(orderSize) ? 0 : negligible + 1
        negligible == 2 && break
    end

    Mmax = isnothing(M) ? something(findlast(>(eps * maximum(orderSize)), orderSize), 1) - 1 : M

    orders = (Nmax + 1 - Mmax):(Nmax + 1 + Mmax)

    A = Aall[orders, :]
    b = ball[orders, :]
    surfaceSize = sizeAll[orders, :]

    # --- the truncation is verified rather than trusted, on the surface of the scatterer: the outgoing radial
    #     functions decrease outward, so that the size of a mode there bounds its contribution everywhere
    tail = maximum(view(surfaceSize, :, max(Nmax, 1):(Nmax + 1))) / maximum(surfaceSize)

    return SpheroidalModes{T}(c, Mmax, Nmax, A, b), tail
end


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
