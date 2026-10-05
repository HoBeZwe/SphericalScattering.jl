
# The supertype of the spheroids, `BC` being the condition on their surface. A spheroid is a surface ξ = ξ₀ of its
# spheroidal coordinates, in which the Helmholtz equation separates, so that the solution is one and the same
# series for every shape: the modal machinery is written once for `Spheroid` and accesses the geometry only
# through the functions defined per shape below, `spheroidalCoordinates`, `cartesianCoordinates`,
# `spheroidalMetric`, `spheroidalBasis`, `shape`, `normalizedCircumradius` and `projectionCoordinate`, together
# with `equatorialRadius`. A spheroid is not to be confused with a `Sphere`, which is not a special case of it:
# the spheroidal coordinates degenerate for a sphere.
#
# The type is documented by its constructor, as a docstring of the type itself would be included twice in the
# documentation. Its conditions are restricted to the acoustic ones so far, lest it overlap the electromagnetic
# scatterers.
abstract type Spheroid{BC<:AcousticBoundary} <: Scatterer{BC} end


"""
    OblateSpheroid{BC,R} <: Spheroid{BC}

Oblate spheroid centered in the origin: the surface ``ξ = ξ_0`` of the oblate spheroidal coordinates with the
semifocal distance `semifocal`, rotationally symmetric about `axis`, on which the condition `boundary` holds.

For ``ξ_0 = 0`` it degenerates into a disc, see [`Disc`](@ref).
"""
struct OblateSpheroid{BC<:AcousticBoundary,R} <: Spheroid{BC}
    semifocal::R
    ξ₀::R
    axis::SVector{3,R}
    boundary::BC
end

"""
    OblateSpheroid{BC}(
        equatorialRadius = error("missing argument `equatorialRadius`"),
        polarRadius      = error("missing argument `polarRadius`"),
        axis             = SVector(0.0, 0.0, 1.0)
    )

Constructor for an oblate spheroid, ``a > b``, see [`Spheroid`](@ref).
"""
function OblateSpheroid{BC}(;
    equatorialRadius=error("missing argument `equatorialRadius`"),
    polarRadius=error("missing argument `polarRadius`"),
    axis=SVector(0.0, 0.0, 1.0),
) where {BC<:AcousticBoundary}

    a, b = promote(equatorialRadius, polarRadius)

    b >= 0 || error("The polar radius must not be negative.")
    a > b || error("The equatorial radius must be larger than the polar radius: an oblate spheroid is required.")

    f = sqrt(a^2 - b^2)

    axisNormalized = normalize(SVector{3}(promote(axis...)))

    return OblateSpheroid{BC,typeof(f)}(f, b / f, axisNormalized, BC())
end


"""
    ProlateSpheroid{BC,R} <: Spheroid{BC}

Prolate spheroid centered in the origin: the surface ``ξ = ξ_0 > 1`` of the prolate spheroidal coordinates with
the semifocal distance `semifocal`, rotationally symmetric about `axis`, on which the condition `boundary` holds.

The degenerate surface ``ξ_0 = 1``, the segment of the axis between the foci, is excluded.
"""
struct ProlateSpheroid{BC<:AcousticBoundary,R} <: Spheroid{BC}
    semifocal::R
    ξ₀::R
    axis::SVector{3,R}
    boundary::BC
end

"""
    ProlateSpheroid{BC}(
        equatorialRadius = error("missing argument `equatorialRadius`"),
        polarRadius      = error("missing argument `polarRadius`"),
        axis             = SVector(0.0, 0.0, 1.0)
    )

Constructor for a prolate spheroid, ``b > a > 0``, see [`Spheroid`](@ref).
"""
function ProlateSpheroid{BC}(;
    equatorialRadius=error("missing argument `equatorialRadius`"),
    polarRadius=error("missing argument `polarRadius`"),
    axis=SVector(0.0, 0.0, 1.0),
) where {BC<:AcousticBoundary}

    a, b = promote(equatorialRadius, polarRadius)

    a > 0 || error("The equatorial radius must be positive: a prolate spheroid degenerating into a line segment is not supported.")
    b > a || error("The polar radius must be larger than the equatorial radius: a prolate spheroid is required.")

    f = sqrt((b - a) * (b + a))

    axisNormalized = normalize(SVector{3}(promote(axis...)))

    return ProlateSpheroid{BC,typeof(f)}(f, b / f, axisNormalized, BC())
end


"""
    Spheroid{BC}(
        equatorialRadius = error("missing argument `equatorialRadius`"),
        polarRadius      = error("missing argument `polarRadius`"),
        axis             = SVector(0.0, 0.0, 1.0)
    )

Constructor for a spheroid centered in the origin, where `BC` is [`SoundHard`](@ref) or [`SoundSoft`](@ref).

The spheroid is specified by its `equatorialRadius` ``a`` and its `polarRadius` ``b``, the latter being measured
along the `axis` of revolution, which is normalized. The shape follows from the radii: for ``a > b`` an
[`OblateSpheroid`](@ref) is returned, with the semifocal distance ``f = \\sqrt{a^2 - b^2}`` and the radial
coordinate ``ξ_0 = b / f`` of the surface; for ``b > a`` a [`ProlateSpheroid`](@ref), with
``f = \\sqrt{b^2 - a^2}`` and again ``ξ_0 = b / f``.

For `polarRadius = 0` the oblate spheroid degenerates into a disc of radius ``a``, see [`Disc`](@ref). The
prolate spheroid degenerating into a line segment, `equatorialRadius = 0`, is not supported.

!!! note
    ``a ≠ b`` is required: the spheroidal coordinates degenerate for a sphere, as ``f → 0`` and ``ξ_0 → ∞`` in
    that limit. Use [`HardSphere`](@ref) or [`SoftSphere`](@ref) for a sphere.
"""
function Spheroid{BC}(;
    equatorialRadius=error("missing argument `equatorialRadius`"),
    polarRadius=error("missing argument `polarRadius`"),
    axis=SVector(0.0, 0.0, 1.0),
) where {BC<:AcousticBoundary}

    equatorialRadius == polarRadius &&
        error("Equal radii describe a sphere, for which the spheroidal coordinates degenerate: use a `HardSphere` or a `SoftSphere`.")

    Shape = equatorialRadius > polarRadius ? OblateSpheroid : ProlateSpheroid

    return Shape{BC}(; equatorialRadius=equatorialRadius, polarRadius=polarRadius, axis=axis)
end


"""
    Disc(BC; radius = error("missing argument `radius`"))

Constructor for a disc of radius ``a``, that is, for the degenerate oblate spheroid with vanishing polar radius,
where `BC` is [`SoundHard`](@ref) or [`SoundSoft`](@ref).

Its two faces are distinguished by the sign of ``η``; the surface is ``ξ_0 = 0``.
"""
Disc(::Type{BC}; radius=error("missing argument `radius`"), axis=SVector(0.0, 0.0, 1.0)) where {BC<:AcousticBoundary} =
    OblateSpheroid{BC}(; equatorialRadius=radius, polarRadius=zero(radius), axis=axis)



# --- the geometry of the oblate spheroid, which is all the modal machinery needs to know about the shape

"""
    spheroidalCoordinates(sphere::Spheroid, point)

Convert a point, given in Cartesian coordinates of the frame of the spheroid, see [`frame`](@ref), into the
spheroidal coordinates ``(ξ, η, φ)`` of its shape.
"""
spheroidalCoordinates(sphere::OblateSpheroid, point) = cart2obl(point, sphere.semifocal)

"""
    cartesianCoordinates(sphere::Spheroid, coordinates)

Convert spheroidal coordinates ``(ξ, η, φ)`` of the shape of the spheroid into Cartesian coordinates of its
frame, see [`frame`](@ref).
"""
cartesianCoordinates(sphere::OblateSpheroid, coordinates) = obl2cart(coordinates, sphere.semifocal)

"""
    spheroidalMetric(sphere::Spheroid, coordinates)

Compute the metric coefficients ``(h_ξ, h_η, h_φ)`` of the spheroidal coordinates at the given ones.
"""
spheroidalMetric(sphere::OblateSpheroid, coordinates) = oblateMetric(coordinates, sphere.semifocal)

"""
    spheroidalBasis(sphere::Spheroid, coordinates)

Compute the unit vectors ``(\\hat{e}_ξ, \\hat{e}_η, \\hat{e}_φ)`` of the spheroidal coordinates at the given ones,
in the frame of the spheroid.
"""
spheroidalBasis(sphere::OblateSpheroid, coordinates) = oblateBasis(coordinates, sphere.semifocal)

"""
    shape(sphere::Spheroid)

Returns the shape of the spheroid as the symbol the spheroidal wave functions expect, `:oblate` or `:prolate`.
"""
shape(sphere::OblateSpheroid) = :oblate

"""
    normalizedCircumradius(sphere::Spheroid)

Returns the radius of the sphere circumscribing the spheroid in units of the semifocal distance, the larger of the
two semi-axes: ``\\sqrt{1 + ξ_0^2}`` for an oblate and ``ξ_0`` for a prolate spheroid. Multiplied by the spheroidal
parameter ``c = kf`` it is the counterpart of ``ka`` for a sphere, which bounds the degree, see [`modes`](@ref).
"""
normalizedCircumradius(sphere::OblateSpheroid) = sqrt(1 + sphere.ξ₀^2)

"""
    projectionCoordinate(sphere::Spheroid)

Returns the default radial coordinate of the surface on which the incident field is projected for an excitation
without a known expansion, see [`projectedCoefficients`](@ref).

For an oblate spheroid it is the surface of the scatterer, or ``ξ = 0.5`` for a scatterer flatter than that: at
the degenerate surface ``ξ = 0`` the regular radial functions of odd ``n - m`` vanish, by which the projection
divides, so that a disc has to be projected off its surface. For a prolate spheroid it is the surface itself,
which never degenerates, as ``ξ_0 > 1``; even a thin one close to the excluded line segment is projected there
without loss of accuracy.
"""
projectionCoordinate(sphere::OblateSpheroid) = max(sphere.ξ₀, oftype(sphere.ξ₀, 0.5))



# --- the geometry of the prolate spheroid

spheroidalCoordinates(sphere::ProlateSpheroid, point) = cart2prol(point, sphere.semifocal)

cartesianCoordinates(sphere::ProlateSpheroid, coordinates) = prol2cart(coordinates, sphere.semifocal)

spheroidalMetric(sphere::ProlateSpheroid, coordinates) = prolateMetric(coordinates, sphere.semifocal)

spheroidalBasis(sphere::ProlateSpheroid, coordinates) = prolateBasis(coordinates, sphere.semifocal)

shape(sphere::ProlateSpheroid) = :prolate

normalizedCircumradius(sphere::ProlateSpheroid) = sphere.ξ₀ # the polar semi-axis is the larger one

projectionCoordinate(sphere::ProlateSpheroid) = sphere.ξ₀



"""
    equatorialRadius(sp::Spheroid)

Returns the equatorial radius of the spheroid, ``f \\sqrt{1 + ξ_0^2}`` for an oblate and ``f \\sqrt{ξ_0^2 - 1}`` for
a prolate one.
"""
equatorialRadius(sp::OblateSpheroid) = sp.semifocal * sqrt(1 + sp.ξ₀^2)

equatorialRadius(sp::ProlateSpheroid) = sp.semifocal * sqrt((sp.ξ₀ - 1) * (sp.ξ₀ + 1))


"""
    polarRadius(sp::Spheroid)

Returns the polar radius ``f ξ_0`` of the spheroid, which vanishes for a disc.
"""
polarRadius(sp::Spheroid) = sp.semifocal * sp.ξ₀


"""
    isdisc(sp::Spheroid)

Returns whether the spheroid is degenerate, that is, a disc.
"""
isdisc(sp::Spheroid) = false

isdisc(sp::OblateSpheroid) = iszero(sp.ξ₀)


"""
    frame(sphere::Spheroid)

Compute the rotation matrix whose third column is the axis of the spheroid.

It maps the spheroid frame, in which the axis of revolution is ``\\hat{e}_z`` and the spheroidal coordinates are
defined, to the global frame. Since the pressure is a scalar, an arbitrarily oriented spheroid requires no more
than transforming the points into that frame; nothing has to be rotated back.
"""
function frame(sphere::Spheroid)

    â = sphere.axis
    R = eltype(â)

    â == SVector{3,R}(0, 0, 1) && return SMatrix{3,3,R}(I)

    # --- right-handed rotation taking ê_z to â via the Rodrigues formula
    aux = SVector{3,R}(0, 0, 1) × â
    rotAxis = iszero(norm(aux)) ? SVector{3,R}(0, 1, 0) : normalize(aux)

    cosϑ = â[3]
    sinϑ = sqrt(â[1]^2 + â[2]^2)

    K = SMatrix{3,3,R}([
        0 -rotAxis[3] rotAxis[2]
        rotAxis[3] 0 -rotAxis[1]
        -rotAxis[2] rotAxis[1] 0
    ])

    return I + sinϑ * K + (1 - cosϑ) * K * K
end


"""
    isinside(scatterer::Spheroid, point)

Returns whether the point lies inside the spheroid, that is, whether its radial coordinate is smaller than that
of the surface, see [`isinside`](@ref).
"""
isinside(scatterer::Spheroid, point) = spheroidalCoordinates(scatterer, frame(scatterer)' * point)[1] < scatterer.ξ₀


"""
    outwardNormal(sphere::Spheroid, point)

Compute the outward unit normal of the spheroid at the surface point having the direction of `point`.

The normal is ``\\hat{e}_ξ`` of the spheroidal coordinates, which is not parallel to the position vector unless
the scatterer is a sphere: for a disc it is ``\\pm \\hat{e}_z``, the sign following the face. The point and the
returned normal are in Cartesian coordinates of the global frame.
"""
function outwardNormal(sphere::Spheroid, point)

    R = frame(sphere)

    ξ, η, φ = spheroidalCoordinates(sphere, R' * point)

    êξ, ~, ~ = spheroidalBasis(sphere, SVector(sphere.ξ₀, η, φ))

    return R * êξ
end


"""
    outwardNormals(sphere::Spheroid, locations)

Compute the outward unit normals of the spheroid at all `locations`, see [`outwardNormal`](@ref).
"""
outwardNormals(sphere::Spheroid, locations) = map(point -> outwardNormal(sphere, point), locations)
