
"""
    translate(points, translation::SVector{3,T})

Translate the points in the direction of the translation vector. 

All inputs are assumed to be in Cartesian coordinates.
"""
function translate(points, translation::SVector{3,T}) where {T}

    translation == SVector{3,T}(0.0, 0.0, 0.0) && return points # no translation

    points_shifted = similar(points)
    for (ind, p) in enumerate(points)
        points_shifted[ind] = p + translation
    end

    return points_shifted
end



"""
    rotate(excitation::Excitation, vectors_list; inverse=false)

Determine rotation matrix and perform rotation for general excitations. 

The points are assumed to be in Cartesian coordinates.

The vectors_list is NOT modified.
"""
function rotate(excitation::Excitation, vectors_list; inverse=false)

    vectors_list_rot = deepcopy(vectors_list)
    rotate!(excitation, vectors_list_rot; inverse=inverse)

    return vectors_list_rot
end


"""
    rotate!(excitation::Excitation, vectors_list; inverse=false)

Determine rotation matrix and perform rotation for a general excitation. 

The points are assumed to be in Cartesian coordinates.

The vectors_list IS modified (overwritten).
"""
function rotate!(excitation::Excitation, vectors_list; inverse=false)

    # --- rotation matrix
    R = rotationMatrix(excitation)
    isnothing(R) && return nothing

    # --- inverse?
    inverse && (R = R')             # for inverse take transpose: inv(R) = R' 

    # --- perform rotation
    for ind in eachindex(vectors_list)
        vectors_list[ind] = R * vectors_list[ind]
    end

    return nothing
end


"""
    rotationMatrix(excitation::PlaneWave)

Determine rotation matrix for a plane wave excitation.
"""
function rotationMatrix(excitation::PlaneWave)

    p = excitation.polarization # unit vector
    k = excitation.direction    # unit vector
    T = eltype(p)

    # --- k == ez and p == ex  => no ratation
    k == SVector{3,T}(0, 0, 1) && p == SVector{3,T}(1, 0, 0) && return nothing

    # --- rotation matrix
    R = SMatrix{3,3,T}([p k × p k])

    return R
end


"""
    rotationMatrix(excitation::Union{Dipole,RingCurrent})

Determine rotation matrix for a dipole or a ring-current excitation via the Rodrigues formula. 
"""
function rotationMatrix(excitation::Union{Dipole,RingCurrent,UniformField})

    p = orientation(excitation) # unit vector
    T = eltype(p)

    # --- p == ez  => no rotation
    p == SVector{3,T}(0, 0, 1) && return nothing

    # --- rotation axis
    aux = SVector{3,T}(0, 0, 1) × p
    if norm(aux) == T(0) # rotation to -ez         
        rotAxis = SVector{3,T}(0, 1, 0) # rotate around y-axis
    else # all other cases
        rotAxis = normalize(aux)
    end

    cosϑ = p[3]                     # inner product of unit vectors:  ez ⋅ p  = p[3] = cos(ϑ)
    sinϑ = sqrt(p[1]^2 + p[2]^2)    # cross product of unit vectors: |ez × p| = sin(ϑ)

    # --- auxiliary matrix
    K = SMatrix{3,3,T}([
         0          -rotAxis[3]  rotAxis[2]
         rotAxis[3]  0          -rotAxis[1]
        -rotAxis[2]  rotAxis[1]  0
    ])

    # --- put it together
    R = I + sinϑ * K + (1 - cosϑ) * K * K  # Rodriguez rotation formula

    return R
end



"""
    convertSpherical2Cartesian(F_sph, point_sph)

Takes a 3D (3 entry) vector `F_sph` and converts it from its spherical basis to its Cartesian basis representation.

The location of the vector has to be provided in spherical coordinates by `point_sph` ordered as ``(r, ϑ, φ)``.
"""
function convertSpherical2Cartesian(F_sph, point_sph)

    T = eltype(F_sph)

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    return SVector{3,T}(
        F_sph[1] * sinϑ * cosϕ + F_sph[2] * cosϑ * cosϕ - F_sph[3] * sinϕ,
        F_sph[1] * sinϑ * sinϕ + F_sph[2] * cosϑ * sinϕ + F_sph[3] * cosϕ,
        F_sph[1] * cosϑ - F_sph[2] * sinϑ,
    )
end


"""
convertSpherical2Cartesian(F_sph, point_sph)

Takes a 3D (3 entry) vector `F_cart` and converts it from its Cartesian basis to its spherical basis representation.

The location of the vector has to be provided in spherical coordinates by `point_sph` ordered as ``(r, ϑ, φ)``.
"""

function convertCartesian2Spherical(F_cart, point_sph)

    T = eltype(F_cart)

    sinϑ = sin(point_sph[2])
    cosϑ = cos(point_sph[2])

    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    return SVector{3,T}(
        F_cart[1] * sinϑ * cosϕ + F_cart[2] * sinϑ * sinϕ + F_cart[3] * cosϑ,
        F_cart[1] * cosϑ * cosϕ + F_cart[2] * cosϑ * sinϕ - F_cart[3] * sinϑ,
        -F_cart[1] * sinϕ + F_cart[2] * cosϕ,
    )
end


"""
    cart2sph(vec)

Convert a 3D (3 entry) point from Cartesian to spherical coordinates.
"""
function cart2sph(vec)

    x = vec[1]
    y = vec[2]
    z = vec[3]

    T = eltype(x)

    r = hypot(x, y, z)
    ϑ = r ≈ T(0) ? T(0) : acos(z / r)   # ∈ [0, π] with ϑ = 0 for x = y = z = 0 case
    φ = atan(y, x)                      # ∈ [-π, π]

    return SVector{3,T}(r, ϑ, φ)
end


"""
    sph2cart(vec)

Convert a 3D (e entry) point from spherical to Cartesian coordinates.
"""
function sph2cart(vec)

    r = vec[1]
    ϑ = vec[2]
    ϕ = vec[3]

    T = typeof(r)

    x = r * sin(ϑ) * cos(ϕ)
    y = r * sin(ϑ) * sin(ϕ)
    z = r * cos(ϑ)

    return SVector{3,T}(x, y, z)
end


"""
    obl2cart(vec, semifocal)

Convert a 3D (3 entry) point from oblate spheroidal to Cartesian coordinates.

The oblate spheroidal coordinates ``(ξ, η, φ)`` with ``ξ ≥ 0`` and ``η ∈ [-1, 1]`` are defined by
``x = f \\sqrt{(1 + ξ^2)(1 - η^2)} \\cos φ``, ``y = f \\sqrt{(1 + ξ^2)(1 - η^2)} \\sin φ`` and ``z = f ξ η``,
where ``f`` denotes the semifocal distance. The surface ``ξ = ξ_0`` is an oblate spheroid with equatorial radius
``f \\sqrt{1 + ξ_0^2}`` and polar radius ``f ξ_0``; the degenerate surface ``ξ = 0`` is the disc of radius ``f``,
its two faces being distinguished by the sign of ``η``.
"""
function obl2cart(vec, semifocal)

    ξ = vec[1]
    η = vec[2]
    φ = vec[3]

    T = typeof(ξ)

    ρ = semifocal * sqrt((1 + ξ^2) * (1 - η^2))

    return SVector{3,T}(ρ * cos(φ), ρ * sin(φ), semifocal * ξ * η)
end


"""
    cart2obl(vec, semifocal)

Convert a 3D (3 entry) point from Cartesian to oblate spheroidal coordinates ``(ξ, η, φ)``.

Inverting the definition in [`obl2cart`](@ref) leads to the quadratic ``u^2 + (1 - A - B) u - B = 0`` for
``u = ξ^2``, where ``A = (x^2 + y^2) / f^2`` and ``B = z^2 / f^2``, of which the non-negative root is taken.

For a point in the plane ``z = 0`` inside the disc the sign of ``η``, that is, the face of the disc, cannot be
recovered; the positive one is returned. This holds for ``z = -0`` as well, whose sign is an artifact of rounding,
e.g., of the rotation into the frame of a scatterer, rather than a choice of the face.
"""
function cart2obl(vec, semifocal)

    x = vec[1]
    y = vec[2]
    z = vec[3]

    T = eltype(x)

    A = (x^2 + y^2) / semifocal^2
    B = z^2 / semifocal^2

    # --- non-negative root of u² + (1 - A - B) u - B = 0
    u = ((A + B - 1) + sqrt((A + B - 1)^2 + 4 * B)) / 2
    u = max(u, T(0.0))

    # --- η² follows either from B = u η² or from A = (1 + u)(1 - η²). The latter cancels for η² → 0,
    #     the former involves no subtraction at all but degenerates as u → 0, that is, on a disc. Hence
    #     the one is taken where the other would cancel
    v = clamp(1 - A / (1 + u), T(0.0), T(1.0))
    v < T(0.5) && u > eps(T) && (v = clamp(B / u, T(0.0), T(1.0)))

    ξ = sqrt(u)
    η = z < 0 ? -sqrt(v) : sqrt(v) # not `copysign`, which would honour the sign of a vanishing z
    φ = atan(y, x)

    return SVector{3,T}(ξ, η, φ)
end


"""
    oblateMetric(vec, semifocal)

Compute the metric coefficients ``(h_ξ, h_η, h_φ)`` of the oblate spheroidal coordinates at the point `vec`,
which is given in oblate spheroidal coordinates.

They read ``h_ξ = f \\sqrt{(ξ^2 + η^2) / (1 + ξ^2)}``, ``h_η = f \\sqrt{(ξ^2 + η^2) / (1 - η^2)}`` and
``h_φ = f \\sqrt{(1 + ξ^2)(1 - η^2)}``. Since the outward normal of the surface ``ξ = ξ_0`` is
``\\hat{e}_ξ = h_ξ^{-1} ∂\\bm{r} / ∂ξ``, the normal derivative is ``∂/∂n = h_ξ^{-1} ∂/∂ξ``.

Note that ``h_ξ = f |η|`` on the disc ``ξ = 0``, so that the normal derivative is singular at its rim
``η = 0``: this is where the edge singularity of a disc enters.
"""
function oblateMetric(vec, semifocal)

    ξ = vec[1]
    η = vec[2]

    T = typeof(ξ)

    # 1 - η² is formed as (1 - η)(1 + η): close to the poles 1 - η is exact, whereas 1 - η² would cancel and
    # deviate from the η at which the angular functions are evaluated
    sη² = (1 - η) * (1 + η)

    hξ = semifocal * sqrt((ξ^2 + η^2) / (1 + ξ^2))
    hη = semifocal * sqrt((ξ^2 + η^2) / sη²)
    hφ = semifocal * sqrt((1 + ξ^2) * sη²)

    return SVector{3,T}(hξ, hη, hφ)
end



"""
    oblateBasis(vec, semifocal)

Compute the unit vectors ``(\\hat{e}_ξ, \\hat{e}_η, \\hat{e}_φ)`` of the oblate spheroidal coordinates at the
point `vec`, which is given in oblate spheroidal coordinates.

They follow from the tangent vectors ``∂\\bm{r} / ∂ξ``, ``∂\\bm{r} / ∂η`` and ``∂\\bm{r} / ∂φ`` of the forward
transform, normalized by the metric coefficients, see [`oblateMetric`](@ref). The first one is the outward normal
of the surface ``ξ = \\mathrm{const}``.

!!! note
    The basis degenerates at the rim of a disc, ``ξ = 0`` and ``η = 0``, where ``\\hat{e}_ξ`` and
    ``\\hat{e}_η`` vanish, and on the axis, ``|η| = 1``, where ``\\hat{e}_φ`` is not determined. These are the
    edge and the poles of the scatterer.
"""
function oblateBasis(vec, semifocal)

    ξ = vec[1]
    η = vec[2]
    φ = vec[3]

    T = typeof(ξ)

    f = semifocal

    sξ = sqrt(1 + ξ^2)
    sη = sqrt((1 - η) * (1 + η)) # see `oblateMetric`

    h = oblateMetric(vec, semifocal)

    êξ = SVector{3,T}(f * ξ * sη / sξ * cos(φ), f * ξ * sη / sξ * sin(φ), f * η) / h[1]
    êη = SVector{3,T}(-f * η * sξ / sη * cos(φ), -f * η * sξ / sη * sin(φ), f * ξ) / h[2]
    êφ = SVector{3,T}(-sin(φ), cos(φ), T(0.0))

    return êξ, êη, êφ
end



"""
    prol2cart(vec, semifocal)

Convert a 3D (3 entry) point from prolate spheroidal to Cartesian coordinates.

The prolate spheroidal coordinates ``(ξ, η, φ)`` with ``ξ ≥ 1`` and ``η ∈ [-1, 1]`` are defined by
``x = f \\sqrt{(ξ^2 - 1)(1 - η^2)} \\cos φ``, ``y = f \\sqrt{(ξ^2 - 1)(1 - η^2)} \\sin φ`` and ``z = f ξ η``,
where ``f`` denotes the semifocal distance. The surface ``ξ = ξ_0`` is a prolate spheroid with equatorial radius
``f \\sqrt{ξ_0^2 - 1}`` and polar radius ``f ξ_0``; the degenerate surface ``ξ = 1`` is the segment of the axis
between the foci.
"""
function prol2cart(vec, semifocal)

    ξ = vec[1]
    η = vec[2]
    φ = vec[3]

    T = typeof(ξ)

    # the differences are formed as products, which do not cancel close to the axis
    ρ = semifocal * sqrt((ξ - 1) * (ξ + 1) * (1 - η) * (1 + η))

    return SVector{3,T}(ρ * cos(φ), ρ * sin(φ), semifocal * ξ * η)
end


"""
    cart2prol(vec, semifocal)

Convert a 3D (3 entry) point from Cartesian to prolate spheroidal coordinates ``(ξ, η, φ)``.

Inverting the definition in [`prol2cart`](@ref) leads to the quadratic ``u^2 - (1 + A + B) u + B = 0`` for
``u = ξ^2``, where ``A = (x^2 + y^2) / f^2`` and ``B = z^2 / f^2``, of which the larger root is taken, as
``u ≥ 1``. Its discriminant is written as ``(1 - B)^2 + A (A + 2 + 2B)``, a sum of non-negative terms.

On the axis, ``x = y = 0``, the sign of ``η`` is that of ``z``; for ``z = 0`` the positive one is returned, as
for [`cart2obl`](@ref).
"""
function cart2prol(vec, semifocal)

    x = vec[1]
    y = vec[2]
    z = vec[3]

    T = eltype(x)

    A = (x^2 + y^2) / semifocal^2
    B = z^2 / semifocal^2

    u = ((1 + A + B) + sqrt((1 - B)^2 + A * (A + 2 + 2 * B))) / 2

    # --- η² follows either from B = u η² or from A = (u - 1)(1 - η²). The former involves no subtraction but is
    #     not exact on the axis, where the latter yields η² = 1 exactly; the latter cancels for η² → 0, though,
    #     and degenerates on the segment between the foci, u = 1. Hence each is taken where the other is worse
    v = clamp(B / u, T(0.0), T(1.0))
    v > T(0.5) && u - 1 > eps(T) && (v = clamp(1 - A / (u - 1), T(0.0), T(1.0)))

    ξ = sqrt(u)
    η = z < 0 ? -sqrt(v) : sqrt(v)
    φ = atan(y, x)

    return SVector{3,T}(ξ, η, φ)
end


"""
    prolateMetric(vec, semifocal)

Compute the metric coefficients ``(h_ξ, h_η, h_φ)`` of the prolate spheroidal coordinates at the point `vec`,
which is given in prolate spheroidal coordinates.

They read ``h_ξ = f \\sqrt{(ξ^2 - η^2) / (ξ^2 - 1)}``, ``h_η = f \\sqrt{(ξ^2 - η^2) / (1 - η^2)}`` and
``h_φ = f \\sqrt{(ξ^2 - 1)(1 - η^2)}``, the outward normal of the surface ``ξ = ξ_0`` being
``\\hat{e}_ξ = h_ξ^{-1} ∂\\bm{r} / ∂ξ``. As for [`oblateMetric`](@ref), the differences are formed as products.
"""
function prolateMetric(vec, semifocal)

    ξ = vec[1]
    η = vec[2]

    T = typeof(ξ)

    sξ² = (ξ - 1) * (ξ + 1)
    sη² = (1 - η) * (1 + η)
    d² = (ξ - η) * (ξ + η)

    hξ = semifocal * sqrt(d² / sξ²)
    hη = semifocal * sqrt(d² / sη²)
    hφ = semifocal * sqrt(sξ² * sη²)

    return SVector{3,T}(hξ, hη, hφ)
end


"""
    prolateBasis(vec, semifocal)

Compute the unit vectors ``(\\hat{e}_ξ, \\hat{e}_η, \\hat{e}_φ)`` of the prolate spheroidal coordinates at the
point `vec`, which is given in prolate spheroidal coordinates, see [`oblateBasis`](@ref).

!!! note
    The basis degenerates on the axis, ``|η| = 1``, where ``\\hat{e}_φ`` is not determined: at the poles of the
    scatterer.
"""
function prolateBasis(vec, semifocal)

    ξ = vec[1]
    η = vec[2]
    φ = vec[3]

    T = typeof(ξ)

    f = semifocal

    sξ = sqrt((ξ - 1) * (ξ + 1))
    sη = sqrt((1 - η) * (1 + η))

    h = prolateMetric(vec, semifocal)

    êξ = SVector{3,T}(f * ξ * sη / sξ * cos(φ), f * ξ * sη / sξ * sin(φ), f * η) / h[1]
    êη = SVector{3,T}(-f * η * sξ / sη * cos(φ), -f * η * sξ / sη * sin(φ), f * ξ) / h[2]
    êφ = SVector{3,T}(-sin(φ), cos(φ), T(0.0))

    return êξ, êη, êφ
end
