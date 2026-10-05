
# [Sound-Hard/Soft Sphere](@id acScattererAPI)

The acoustic counterparts of the PEC sphere have radius ``r`` and are assumed to be located in the origin. On a sound-hard surface the normal velocity, and hence the normal derivative of the total pressure, vanishes; on a sound-soft (pressure release) surface the total pressure vanishes.

## API

```@docs
HardSphere
SoftSphere
```

Both are aliases: a sound-hard sphere is a `Sphere{SoundHard}`, that is, a sphere carrying the boundary condition [`SoundHard`](@ref) as its type parameter. Every scatterer is a [`Scatterer`](@ref) parametrized by the condition on its surface, which allows to address, e.g., all acoustic scatterers at once as `Scatterer{<:AcousticBoundary}`.

!!! note
    These two are the only acoustic scatterers with a closed surface of revolution for which the cheap spherical series applies. The spheroids and the disc are documented separately, together with the accuracy of their solution, under [Spheroid and Disc](@ref ACspheroidAPI).
