
"""
    AcousticBoundary

Supertype of the acoustic boundary conditions, which are orthogonal to the geometry of the scatterer and are,
therefore, carried as a type parameter.
"""
abstract type AcousticBoundary end

"""
    SoundHard

The normal velocity, and hence the normal derivative of the total pressure, vanishes on the surface.
"""
struct SoundHard <: AcousticBoundary end

"""
    SoundSoft

The total pressure vanishes on the surface (pressure release).
"""
struct SoundSoft <: AcousticBoundary end



struct Spheroid{BC<:AcousticBoundary,R} <: Sphere
    semifocal::R
    ξ₀::R
    axis::SVector{3,R}
end

"""
    Spheroid{BC}(
        equatorialRadius = error("missing argument `equatorialRadius`"),
        polarRadius      = error("missing argument `polarRadius`"),
        axis             = SVector(0.0, 0.0, 1.0)
    )

Constructor for an oblate spheroid centered in the origin, where `BC` is [`SoundHard`](@ref) or
[`SoundSoft`](@ref).

The spheroid is specified by its `equatorialRadius` ``a`` and its `polarRadius` ``b``, from which the semifocal
distance ``f = \\sqrt{a^2 - b^2}`` and the radial coordinate ``ξ_0 = b / f`` of the surface are determined. The
`axis` is the axis of revolution and is normalized.

For `polarRadius = 0` the spheroid degenerates into a disc of radius ``a``, see [`Disc`](@ref).

!!! note
    Strictly ``a > b`` is required: the oblate spheroidal coordinates degenerate for a sphere, as ``f → 0`` and
    ``ξ_0 → ∞`` in that limit. Use [`HardSphere`](@ref) or [`SoftSphere`](@ref) for a sphere.
"""
function Spheroid{BC}(;
    equatorialRadius=error("missing argument `equatorialRadius`"),
    polarRadius=error("missing argument `polarRadius`"),
    axis=SVector(0.0, 0.0, 1.0),
) where {BC<:AcousticBoundary}

    a, b = promote(equatorialRadius, polarRadius)

    b >= 0 || error("The polar radius must not be negative.")
    a > b || error("The equatorial radius must be larger than the polar radius: an oblate spheroid is required.")

    f = sqrt(a^2 - b^2)

    axisNormalized = normalize(SVector{3}(promote(axis...)))

    return Spheroid{BC,typeof(f)}(f, b / f, axisNormalized)
end


"""
    Disc(BC; radius = error("missing argument `radius`"))

Constructor for a disc of radius ``a``, that is, for the degenerate oblate spheroid with vanishing polar radius,
where `BC` is [`SoundHard`](@ref) or [`SoundSoft`](@ref).

Its two faces are distinguished by the sign of ``η``; the surface is ``ξ_0 = 0``.
"""
Disc(::Type{BC}; radius=error("missing argument `radius`"), axis=SVector(0.0, 0.0, 1.0)) where {BC<:AcousticBoundary} =
    Spheroid{BC}(; equatorialRadius=radius, polarRadius=zero(radius), axis=axis)


"""
    equatorialRadius(sp::Spheroid)

Returns the equatorial radius ``f \\sqrt{1 + ξ_0^2}`` of the spheroid.
"""
equatorialRadius(sp::Spheroid) = sp.semifocal * sqrt(1 + sp.ξ₀^2)


"""
    polarRadius(sp::Spheroid)

Returns the polar radius ``f ξ_0`` of the spheroid, which vanishes for a disc.
"""
polarRadius(sp::Spheroid) = sp.semifocal * sp.ξ₀


"""
    isdisc(sp::Spheroid)

Returns whether the spheroid is degenerate, that is, a disc.
"""
isdisc(sp::Spheroid) = iszero(sp.ξ₀)
"""
    frame(sphere::Spheroid)

Compute the rotation matrix whose third column is the axis of the spheroid.

It maps the spheroid frame, in which the axis of revolution is ``\\hat{e}_z`` and the oblate spheroidal
coordinates are defined, to the global frame. Since the pressure is a scalar, an arbitrarily oriented spheroid
requires no more than transforming the points into that frame; nothing has to be rotated back.
"""
function frame(sphere::Spheroid{BC,R}) where {BC,R}

    â = sphere.axis

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
isinside(scatterer::Spheroid, point) = cart2obl(frame(scatterer)' * point, scatterer.semifocal)[1] < scatterer.ξ₀
