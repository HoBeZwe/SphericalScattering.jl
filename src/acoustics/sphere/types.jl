
struct HardSphere{R} <: Sphere
    radius::R
end

"""
    HardSphere(
        radius = error("missing argument `radius`")
    )

Constructor for a sound-hard sphere.
"""
HardSphere(; radius=error("missing argument `radius`")) = HardSphere(radius)



struct SoftSphere{R} <: Sphere
    radius::R
end

"""
    SoftSphere(
        radius = error("missing argument `radius`")
    )

Constructor for a sound-soft sphere.
"""
SoftSphere(; radius=error("missing argument `radius`")) = SoftSphere(radius)


"""
    isinside(scatterer::Union{HardSphere,SoftSphere}, point)

Returns whether the point lies inside the sphere, see [`isinside`](@ref).
"""
isinside(scatterer::Union{HardSphere,SoftSphere}, point) = norm(point) < scatterer.radius
