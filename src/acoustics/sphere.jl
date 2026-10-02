
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