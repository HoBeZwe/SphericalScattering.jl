
"""
    PEC

Perfect electric conductor: the tangential electric field vanishes on the surface.
"""
struct PEC <: ElectromagneticBoundary end


"""
    Dielectric{C} <: ElectromagneticBoundary

The surface encloses a homogeneous, isotropic `filling`, see [`Medium`](@ref): the tangential electric and
magnetic fields are continuous across it.
"""
struct Dielectric{C} <: ElectromagneticBoundary
    filling::Medium{C}
end


"""
    ThinImpedanceLayer{R,C} <: ElectromagneticBoundary

The surface encloses a homogeneous, isotropic `filling` coated by a thin layer of the medium `thinlayer` and of
the given `thickness`, see [`Medium`](@ref).

The condition is an approximate one: the displacement field in the layer is assumed to be radial, see
[`DielectricSphereThinImpedanceLayer`](@ref).
"""
struct ThinImpedanceLayer{R,C} <: ElectromagneticBoundary
    thickness::R
    thinlayer::Medium{C}
    filling::Medium{C}
end


"""
    Layered{Core,N,R,C} <: ElectromagneticBoundary

The surface encloses `N` concentric shells around a core, on whose surface the condition `core` holds, a
[`Dielectric`](@ref) or a [`PEC`](@ref) one.

Seen from outside, the shells form a condition on the outer surface, whose radius is that of the
[`Sphere`](@ref) carrying it. Hence only the inner interfaces are stored, ascendingly: the shell `k` extends
from `radii[k]` to `radii[k + 1]`, the outermost one to the radius of the sphere, and is filled by
`filling[k]`. The core fills the interior of `radii[1]`, or that of the sphere if there are no shells.
"""
struct Layered{Core<:ElectromagneticBoundary,N,R,C} <: ElectromagneticBoundary
    radii::SVector{N,R}
    filling::SVector{N,Medium{C}}
    core::Core
end
