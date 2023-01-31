
"""
    translate(points, translation::SVector{3,T})

Translate the points in negative! direction of the translation vector. 

All inputs are assumed do be in Cartesian coordinates.
"""
function translate(points, translation::SVector{3,T}) where {T}

    translation == SVector{3,T}(0.0, 0.0, 0.0) && return points # no translation

    points_shifted = similar(points)
    for (ind, p) in enumerate(points)
        points_shifted[ind] = p - translation
    end

    return points_shifted
end


"""
    rotate!(points, rotation::SVector{2,T})

Rotate points around x-axis and y-axis away from the z-axis by angles (in radian) specified in rotation[1] and rotation[2]. 

The points are assumed do be in Cartesian coordinates.
"""
function rotate!(points, rotation::SVector{2,T}) where {T}

    rotation == SVector{2,T}(0.0, 0.0) && return nothing # no rotation required

    # --- perform rotation
    # wx = rotation[1] # rotation angle around x-axis (away from z-axis)
    # wy = rotation[2] # rotation angle around y-axis (away from z-axis)

    # cosy = cos(wy)
    # siny = sin(wy)
    # cosx = cos(wx)
    # sinx = sin(wx)

    # Ry = [
    #     cosy 0 siny
    #     0 1 0
    #     -siny 0 cosy
    # ]

    # Rx = [
    #     1 0 0
    #     0 cosx sinx
    #     0 -sinx cosx
    # ]

    # R = Rx * Ry

    # for (ind, p) in enumerate(points)
    #     points[ind] = R * p
    # end

    return nothing
end



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


# ----------------------- Convert the computed field from Cartesian to spherical representation
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


function cart2sph(vec)
    x = vec[1]
    y = vec[2]
    z = vec[3]

    T = eltype(x)

    r = sqrt(x^2 + y^2 + z^2)
    if x ≈ T(0) && y ≈ T(0)
        if z >= 0
            t = T(0)
        else
            t = T(π)
        end
    else
        t = (atan(sqrt(x^2 + y^2), z) + π) % π  # map to range [0, π]
    end
    p = (atan(y, x) + 2 * π) % 2π               # map to range [0, 2π]

    return SVector{3,T}(r, t, p)
end


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
