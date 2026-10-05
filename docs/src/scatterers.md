
# [Scatterers and Boundary Conditions](@id scatterersConcept)

A scatterer is the combination of a geometry and a condition on its surface. The geometry is the type, the condition its type parameter: a PEC sphere is a `Sphere{PEC}`, a sound-hard oblate spheroid an `OblateSpheroid{SoundHard}`. All of them are subtypes of [`Scatterer{BC}`](@ref Scatterer), where `BC` is the boundary condition.

The two are independent of each other, so that either can be addressed alone:

| addressed by            | example                            | holds for                                  |
|:----------------------- |:---------------------------------- |:------------------------------------------ |
| geometry                | `Sphere`                           | spheres with any condition                 |
| geometry and condition  | `Sphere{SoundHard}`                | sound-hard spheres                         |
| physics                 | `Scatterer{<:AcousticBoundary}`    | all acoustic scatterers                    |


---
## Type Hierarchy

The scatterers:
```
Scatterer{BC}
├── Sphere{BC}               centered in the origin
└── Spheroid{BC}             acoustic conditions only (so far)
    ├── OblateSpheroid{BC}   including the disc
    └── ProlateSpheroid{BC}
```

The boundary conditions, split by physics:
```
Boundary
├── ElectromagneticBoundary
│   ├── PEC
│   ├── Dielectric
│   ├── Layered
│   └── ThinImpedanceLayer
└── AcousticBoundary
    ├── SoundHard
    └── SoundSoft
```

See [`Sphere`](@ref SphericalScattering.Sphere), [`Spheroid`](@ref), [`OblateSpheroid`](@ref), [`ProlateSpheroid`](@ref), [`Boundary`](@ref SphericalScattering.Boundary), [`ElectromagneticBoundary`](@ref), [`AcousticBoundary`](@ref), [`PEC`](@ref), [`Dielectric`](@ref), [`Layered`](@ref), [`ThinImpedanceLayer`](@ref), [`SoundHard`](@ref), and [`SoundSoft`](@ref).

A condition which requires data, such as the filling of a [`Dielectric`](@ref), is stored as a value besides being the type parameter, `sp.boundary`, whereas a parameter-free one, such as [`PEC`](@ref) or [`SoundHard`](@ref), occupies no memory.

!!! note
    A sphere is not a special case of a spheroid: the spheroidal coordinates degenerate for a sphere, see the [oblate](@ref oblateCoords) and the [prolate](@ref prolateCoords) spheroidal coordinates. A disc, in contrast, is the oblate spheroid of vanishing polar radius.


---
## Constructors

Every scatterer is defined by a constructor with keyword arguments:

| scatterer | constructor | type | details |
|:--------- |:----------- |:---- |:------- |
| PEC sphere | `PECSphere(; radius)` | `Sphere{PEC}` | [PEC sphere](@ref pecAPI) |
| dielectric sphere | `DielectricSphere(; radius, filling)` | `Sphere{<:Dielectric}` | [dielectric sphere](@ref dielecAPI) |
| layered dielectric sphere | `LayeredSphere(; radii, filling)` | `Sphere{<:Layered{<:Dielectric}}` | [layered sphere](@ref mlDielecAPI) |
| layered dielectric sphere with PEC core | `LayeredSpherePEC(; radii, filling)` | `Sphere{<:Layered{PEC}}` | [layered sphere with PEC core](@ref mlDielecPecAPI) |
| dielectric sphere with thin impedance layer | `DielectricSphereThinImpedanceLayer(; radius, thickness, thinlayer, filling)` | `Sphere{<:ThinImpedanceLayer}` | [thin impedance layer](@ref dielecimpedAPI) |
| sound-hard sphere | `HardSphere(; radius)` | `Sphere{SoundHard}` | [sound-hard/soft sphere](@ref acScattererAPI) |
| sound-soft sphere | `SoftSphere(; radius)` | `Sphere{SoundSoft}` | [sound-hard/soft sphere](@ref acScattererAPI) |
| spheroid | `Spheroid{BC}(; equatorialRadius, polarRadius, axis)` | `OblateSpheroid{BC}` or `ProlateSpheroid{BC}` | [spheroid and disc](@ref ACspheroidAPI) |
| disc | `Disc(BC; radius, axis)` | `OblateSpheroid{BC}` | [spheroid and disc](@ref ACspheroidAPI) |

The names of the spheres are aliases of the respective `Sphere{BC}`, so that they can be dispatched on as well. A [`Spheroid`](@ref) returns the oblate or the prolate shape, depending on which of its radii is the larger one.

Alternatively, a sphere can be defined by its condition directly, which yields the same scatterers:
```@example scatterers
using SphericalScattering

sp = SphericalScattering.Sphere{SoundHard}(; radius=1.0)   # parameter-free condition: the type suffices

sp == HardSphere(; radius=1.0)
```
```@example scatterers
sp = SphericalScattering.Sphere(; radius=1.0, boundary=Dielectric(Medium(2.0, 1.0)))   # condition requiring data

typeof(sp)
```


---
## Combining Scatterers and Excitations

The physics of the scatterer has to match that of the excitation: an electromagnetic excitation requires a `Scatterer{<:ElectromagneticBoundary}`, an acoustic one a `Scatterer{<:AcousticBoundary}`. Otherwise a descriptive error is thrown:
```@example scatterers
using StaticArrays

ex = planeWave(; frequency=1e8)   # electromagnetic
sp = HardSphere(; radius=1.0)     # acoustic

try
    scatteredfield(sp, ex, ElectricField([SVector(2.0, 0.0, 0.0)]))
catch err
    showerror(stdout, err)
end
```

Which combinations of scatterers and excitations are implemented is summarized in the [Feature Overview](@ref).

The type parameter can also be used to write methods for all scatterers of one physics:
```@example scatterers
physics(sp::Scatterer{<:AcousticBoundary}) = "acoustic"
physics(sp::Scatterer{<:ElectromagneticBoundary}) = "electromagnetic"

physics.([PECSphere(; radius=1.0), Disc(SoundSoft; radius=1.0)])
```


---
## Geometric Queries

Besides the solutions, the following functions answer questions about the geometry:

| function | scatterers | returns |
|:-------- |:---------- |:------- |
| [`isinside(sp, point)`](@ref SphericalScattering.isinside) | spheres and spheroids | whether the point lies inside the scatterer |
| [`equatorialRadius(sp)`](@ref equatorialRadius), [`polarRadius(sp)`](@ref polarRadius) | spheroids | the radii of the spheroid |
| [`isdisc(sp)`](@ref isdisc) | spheroids | whether the spheroid is a disc |
| [`outwardNormal(sp, point)`](@ref outwardNormal), [`outwardNormals(sp, points)`](@ref outwardNormals) | spheroids | the outward normal at the given surface points |
