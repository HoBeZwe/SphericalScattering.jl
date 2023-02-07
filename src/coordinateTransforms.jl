
"""
    translate(points, translation::SVector{3,T})

Translate the points in the direction of the translation vector. 

All inputs are assumed do be in Cartesian coordinates.
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

The points are assumed do be in Cartesian coordinates.

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

The points are assumed do be in Cartesian coordinates.

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
