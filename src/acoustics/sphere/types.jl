
# `HardSphere` and `SoftSphere` are aliases, so that they can be dispatched on as well. Their docstrings are
# attached to the constructors rather than to the aliases, as for the other scatterers: a docstring of the
# binding itself would be included twice in the documentation, by its `@docs` block and by the API reference

const HardSphere = Sphere{SoundHard}

"""
    HardSphere(
        radius = error("missing argument `radius`")
    )

Constructor for a sound-hard sphere, that is, for a [`Sphere`](@ref) with the boundary condition
[`SoundHard`](@ref).

`HardSphere` is an alias of `Sphere{SoundHard}`, so that it can be dispatched on as well.
"""
HardSphere(; radius=error("missing argument `radius`")) = Sphere(radius, SoundHard())


const SoftSphere = Sphere{SoundSoft}

"""
    SoftSphere(
        radius = error("missing argument `radius`")
    )

Constructor for a sound-soft sphere, that is, for a [`Sphere`](@ref) with the boundary condition
[`SoundSoft`](@ref).

`SoftSphere` is an alias of `Sphere{SoundSoft}`, so that it can be dispatched on as well.
"""
SoftSphere(; radius=error("missing argument `radius`")) = Sphere(radius, SoundSoft())
