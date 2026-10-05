
"""
    PressureNormalGradient(sphere::Spheroid, locations)

Construct the Neumann trace for the outward normals of the spheroid at `locations`.

The normals of a spheroid are not parallel to the position vectors, so that the default of
[`PressureNormalGradient`](@ref), the radial direction, is not the outward normal unless the scatterer is a
sphere. This constructor fills in the correct ones, see [`outwardNormal`](@ref).
"""
PressureNormalGradient(sphere::Spheroid, locations) = PressureNormalGradient(locations, outwardNormals(sphere, locations))


"""
    surfaceseries(sphere::Spheroid, md::SpheroidalModes, point, coefficients)

Evaluate the series of the scattered field and its three partial derivatives on the surface of the spheroid, at
the point having the direction of `point`.

Returned is a named tuple `(value, gradient)`, the gradient being given in Cartesian coordinates of the global
frame. It follows from the partial derivatives via
``∇p = h_ξ^{-1} ∂_ξ p \\, \\hat{e}_ξ + h_η^{-1} ∂_η p \\, \\hat{e}_η + h_φ^{-1} ∂_φ p \\, \\hat{e}_φ``, where the
derivative with respect to ``φ`` amounts to a factor ``\\mathrm{j}m``.

Only the direction of `point` is taken into account, the radial coordinate being replaced by that of the surface.

On the axis, ``η = ±1``, the basis vectors ``\\hat{e}_η`` and ``\\hat{e}_φ`` are not determined; the tangential
part of the gradient is obtained as a limit there instead, see [`axialGradient`](@ref).
"""
function surfaceseries(sphere::Spheroid, md::SpheroidalModes{T}, point, coefficients) where {T}

    R = frame(sphere)

    ~, η, φ = spheroidalCoordinates(sphere, R' * point)

    ξ₀ = sphere.ξ₀

    onAxis = isone(abs(η))

    u = Complex{T}(0.0)   # the value and the three partial derivatives
    uξ = Complex{T}(0.0)
    uη = Complex{T}(0.0)
    uφ = Complex{T}(0.0)

    for m in (-md.M):(md.M)

        mAbs = abs(m)

        S = smn(mAbs, mAbs:(md.N), md.c, [η]; spheroid=shape(sphere), normalize=false)
        Rr = rmn(mAbs, mAbs:(md.N), md.c, [ξ₀]; spheroid=shape(sphere), kind=outgoingKind)

        phase = cis(m * φ)

        for (k, n) in enumerate(mAbs:(md.N))
            C = coefficients[m + md.M + 1, n + 1] * phase

            u += C * Rr.value[1, k] * S.value[1, k]
            uξ += C * Rr.derivative[1, k] * S.value[1, k]

            onAxis && continue # the angular derivatives of order one diverge there, see `axialGradient`

            uη += C * Rr.value[1, k] * S.derivative[1, k]
            uφ += C * Rr.value[1, k] * S.value[1, k] * im * m
        end
    end

    h = spheroidalMetric(sphere, SVector(ξ₀, η, φ))
    êξ, êη, êφ = spheroidalBasis(sphere, SVector(ξ₀, η, φ))

    tangential = onAxis ? axialGradient(sphere, md, η, coefficients) : uη / h[2] * êη + uφ / h[3] * êφ

    gradient = R * (uξ / h[1] * êξ + tangential)

    return (value=u, gradient=gradient)
end


"""
    axialGradient(sphere::Spheroid, md::SpheroidalModes, η, coefficients)

Compute the tangential part of the gradient of the series on the axis of the spheroid, ``η = ±1``, in Cartesian
coordinates of its frame.

There ``h_η`` diverges and ``h_φ`` vanishes, while the basis vectors ``\\hat{e}_η`` and ``\\hat{e}_φ`` are not
determined. The tangential gradient is nevertheless finite, as the surface is smooth, and is carried by the
orders ``m = ±1`` alone: their angular functions vanish like ``\\sqrt{1 - η^2}``, those of ``|m| ≥ 2`` faster, and
the order zero has no tangential gradient on the axis. Writing ``S_{1n} = \\sqrt{1 - η^2} \\, g_n``, the two
terms of order ``m`` combine into

```math
\\cfrac{1}{a} \\, C_{mn} \\, R_{1n}(c, ξ_0) \\, g_n(±1) \\, (1, \\mathrm{j}m, 0) \\,,
```

where ``a`` denotes the equatorial radius: the dependence on ``φ`` cancels, as it has to on the axis. The limits ``g_n(±1)`` are extrapolated from two evaluations close to the pole (Richardson), the error being
of the order of the square of their distance from it.
"""
function axialGradient(sphere::Spheroid, md::SpheroidalModes{T}, η, coefficients) where {T}

    G = zero(SVector{3,Complex{T}})

    (md.M < 1 || md.N < 1) && return G # only the order zero is excited, e.g., by an axial plane wave

    δ = T(2)^-26 # 1 - η and 1 + η are exact at the points of the extrapolation

    S(d) = vec(smn(1, 1:(md.N), md.c, [sign(η) * (1 - d)]; spheroid=shape(sphere), normalize=false).value)
    g(d) = S(d) ./ sqrt(d * (2 - d))

    gₙ = 2 .* g(δ) .- g(2δ) # the error of first order in δ cancels

    Rr = rmn(1, 1:(md.N), md.c, [sphere.ξ₀]; spheroid=shape(sphere), kind=outgoingKind).value

    for (k, n) in enumerate(1:(md.N))
        Cp = coefficients[md.M + 2, n + 1] # m = +1
        Cm = coefficients[md.M, n + 1]     # m = -1

        G += Rr[1, k] * gₙ[k] * (Cp * SVector{3,Complex{T}}(1, im, 0) + Cm * SVector{3,Complex{T}}(1, -im, 0))
    end

    return G / equatorialRadius(sphere)
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, quantity::PressureTrace; parameter::Parameter=Parameter())

Compute the Dirichlet trace of the pressure scattered by an oblate spheroid.

It is the scattered pressure evaluated on the surface, that is, at ``ξ = ξ_0``. Only the direction of the point
is taken into account, so that the locations may also be given by the points of a faceted surface mesh.

!!! note
    For a disc the two faces ``η > 0`` and ``η < 0`` share their Cartesian coordinates, so that the face cannot
    be recovered from a location: the one with ``η > 0`` is returned. The other follows from the parity of the
    angular functions, which the degenerate surface selects by the boundary condition:

    | disc | ``γ_0`` | ``γ_1`` |
    |:---- |:------- |:------- |
    | sound-soft | same on both faces | opposite |
    | sound-hard | opposite | vanishes |

    The parity enters because only the modes with even ``n - m`` contribute for a sound-soft disc and only
    those with odd ``n - m`` for a sound-hard one, the angular functions being even and odd in ``η``
    accordingly.
"""
function scatteredfield(
    sphere::Spheroid,
    excitation::AcousticExcitation,
    md::SpheroidalModes,
    point,
    quantity::PressureTrace;
    parameter::Parameter=Parameter(),
)

    return surfaceseries(sphere, md, point, md.A .* md.b).value
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace ``\\hat{n} ⋅ ∇ p_\\mathrm{s}`` of the pressure scattered by an oblate spheroid for the
given `normal`.

In contrast to a sphere, the gradient does not reduce to two terms for an arbitrary normal: the scattered field
of a spheroid depends on ``φ`` as well, so that all three components contribute. For the outward normal
``\\hat{n} = \\hat{e}_ξ`` the tangential ones drop out and the radial derivative remains, see
[`outwardNormal`](@ref) and [`surfaceseries`](@ref).
"""
function scatteredfield(
    sphere::Spheroid,
    excitation::AcousticExcitation,
    md::SpheroidalModes,
    point,
    normal,
    quantity::PressureNormalGradient;
    parameter::Parameter=Parameter(),
)

    return dot(normal, surfaceseries(sphere, md, point, md.A .* md.b).gradient)
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Trace; parameter::Parameter=Parameter())

Compute a surface trace of the pressure scattered by an oblate spheroid at all locations of `quantity`.
"""
function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::PressureTrace; parameter::Parameter=Parameter()
)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    coefficients = md.A .* md.b

    p = progress(length(quantity.locations))

    @tasks for ind in eachindex(quantity.locations)
        F[ind] = surfaceseries(sphere, md, quantity.locations[ind], coefficients).value
        next!(p)
    end
    finish!(p)

    return F
end

function scatteredfield(
    sphere::Spheroid,
    excitation::AcousticExcitation,
    md::SpheroidalModes,
    quantity::PressureNormalGradient;
    parameter::Parameter=Parameter(),
)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    coefficients = md.A .* md.b

    p = progress(length(quantity.locations))

    @tasks for ind in eachindex(quantity.locations)
        F[ind] = dot(quantity.normals[ind], surfaceseries(sphere, md, quantity.locations[ind], coefficients).gradient)
        next!(p)
    end
    finish!(p)

    return F
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::Trace; parameter::Parameter=Parameter())

Compute a surface trace of the pressure scattered by an oblate spheroid, determining the modal coefficients
automatically, see [`modes`](@ref).
"""
function scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::PressureTrace; parameter::Parameter=Parameter())

    md = modes(sphere, excitation; parameter=parameter)

    return scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
end

function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, quantity::PressureNormalGradient; parameter::Parameter=Parameter()
)

    md = modes(sphere, excitation; parameter=parameter)

    return scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
end


"""
    field(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Trace; parameter::Parameter=Parameter())

Compute a total surface trace in the presence of an oblate spheroid, for an incident acoustic excitation.
"""
function field(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::Trace; parameter::Parameter=Parameter()
)

    F = field(excitation, quantity; parameter=parameter)
    F .+= scatteredfield(sphere, excitation, md, quantity; parameter=parameter)

    return F
end


"""
    field(excitation::AcousticExcitation, quantity::PressureJump; parameter::Parameter=Parameter())

Compute the jump of the pressure of an incident acoustic excitation across an open surface, which vanishes.

The incident field is regular across the surface, which is fictitious as far as it is concerned, so that its
jump is zero. Consequently the jump of the total field equals that of the scattered field.
"""
function field(excitation::AcousticExcitation, quantity::PressureJump; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    return zeros(Complex{T}, size(quantity.locations))
end


"""
    checkDisc(sphere::Spheroid)

Ensure that the scatterer is a disc, across which alone a jump is defined.

The check belongs before the loop over the locations: an error thrown inside the parallel loop is wrapped in a
`TaskFailedException` as soon as more than one thread is available, which would make the failure depend on the
number of threads.
"""
checkDisc(sphere::Spheroid) =
    isdisc(sphere) || error("The jump of the pressure is defined across an open surface, that is, across a disc.")


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, point, quantity::PressureJump; parameter::Parameter=Parameter())

Compute the jump ``[p] = p|_+ - p|_-`` of the pressure across a disc.

The two faces of the degenerate surface ``ξ = 0`` are ``η > 0`` and ``η < 0``, and the angular functions obey
``S_{mn}(c, -η) = (-1)^{n-m} S_{mn}(c, η)``, so that

```math
[p] = 2 \\sum_{n-m \\ \\mathrm{odd}} A_{mn} b_{mn} R^{(\\mathrm{out})}_{mn}(c, 0) S_{mn}(c, η)
      \\mathrm{e}^{\\mathrm{j}mφ} \\,,
```

the modes of even ``n - m`` cancelling. Weighting the terms by their parity avoids having to evaluate the two
faces separately, which the Cartesian coordinates of a location cannot distinguish.

Since the degenerate surface retains the modes of odd ``n - m`` for a sound-hard disc and those of even
``n - m`` for a sound-soft one, the jump of the pressure is carried entirely by the sound-hard case and vanishes
identically for the sound-soft one, whose unknown is the jump of the normal derivative instead.

The jump of the total field equals this one, the incident field being continuous across the disc.
"""
function scatteredfield(
    sphere::Spheroid,
    excitation::AcousticExcitation,
    md::SpheroidalModes{T},
    point,
    quantity::PressureJump;
    parameter::Parameter=Parameter(),
) where {T}

    checkDisc(sphere)

    coefficients = md.A .* md.b

    ~, η, φ = spheroidalCoordinates(sphere, frame(sphere)' * point)

    u = Complex{T}(0.0)

    for m in (-md.M):(md.M)

        mAbs = abs(m)

        S = smn(mAbs, mAbs:(md.N), md.c, [η]; spheroid=shape(sphere), normalize=false).value
        Rr = rmn(mAbs, mAbs:(md.N), md.c, [sphere.ξ₀]; spheroid=shape(sphere), kind=outgoingKind).value

        phase = cis(m * φ)

        for (k, n) in enumerate(mAbs:(md.N))
            isodd(n - mAbs) || continue   # the modes of even parity are the same on both faces

            u += 2 * coefficients[m + md.M + 1, n + 1] * Rr[1, k] * S[1, k] * phase
        end
    end

    return u
end


"""
    scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::PressureJump; parameter::Parameter=Parameter())

Compute the jump of the scattered pressure across a disc at all locations of `quantity`.
"""
function scatteredfield(
    sphere::Spheroid, excitation::AcousticExcitation, md::SpheroidalModes, quantity::PressureJump; parameter::Parameter=Parameter()
)

    checkDisc(sphere)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    p = progress(length(quantity.locations))

    @tasks for ind in eachindex(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, md, quantity.locations[ind], quantity; parameter=parameter)
        next!(p)
    end
    finish!(p)

    return F
end


function scatteredfield(sphere::Spheroid, excitation::AcousticExcitation, quantity::PressureJump; parameter::Parameter=Parameter())

    md = modes(sphere, excitation; parameter=parameter)

    return scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
end


"""
    scatteredfield(sphere::Sphere{<:AcousticBoundary}, excitation::AcousticExcitation, quantity::PressureJump; parameter::Parameter=Parameter())

Descriptive error for the jump of the pressure across a closed surface.
"""
function scatteredfield(
    sphere::Sphere{<:AcousticBoundary}, excitation::AcousticExcitation, quantity::PressureJump; parameter::Parameter=Parameter()
)

    return error("The jump of the pressure is defined across an open surface, that is, across a disc.")
end
