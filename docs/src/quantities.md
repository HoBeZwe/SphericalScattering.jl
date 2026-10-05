
# [Quantities](@id quantitiesConcept)

What is computed is selected by a quantity, which wraps the locations of observation:
```julia
F = field(ex, quantity)                # incident field
F = scatteredfield(sp, ex, quantity)   # scattered field
F = field(sp, ex, quantity)            # total field
```
The locations are Cartesian points, e.g., `SVector`s, in an array of any shape, see [Defining Observation Points](@ref). The result has the same shape: per location, a vector quantity yields an `SVector{3}` of its complex Cartesian components, a scalar quantity a complex number.


---
## Overview

| quantity | physics | kind | value | available for |
|:-------- |:------- |:---- |:----- |:------------- |
| `ElectricField` | electromagnetic | field | vector | all electromagnetic excitations |
| `MagneticField` | electromagnetic | field | vector | the time-harmonic excitations |
| `FarField` | both | far field | vector (electromagnetic), scalar (acoustic) | the time-harmonic excitations |
| `ScalarPotential` | electromagnetic | field | scalar | the [uniform static field](@ref uniformEx) |
| `DisplacementField` | electromagnetic | field | vector | the [dielectric sphere with thin impedance layer](@ref dielecimped) |
| `ScalarPotentialJump` | electromagnetic | field | scalar | the [dielectric sphere with thin impedance layer](@ref dielecimped) |
| `Pressure` | acoustic | field | scalar | all acoustic excitations |
| `PressureTrace` | acoustic | trace | scalar | all acoustic excitations |
| [`PressureNormalGradient`](@ref) | acoustic | trace | scalar | all acoustic excitations |
| [`PressureJump`](@ref) | acoustic | trace | scalar | the [disc](@ref ACspheroidAPI) |

Which scatterers are available for which excitation is summarized in the [Feature Overview](@ref). The [radar cross section](@ref rcsPW) is derived from the far field and has a function of its own.


---
## Fields

The fields are evaluated at the given locations. Inside a scatterer, the total field vanishes where the scatterer is impenetrable, that is, inside a PEC, a sound-hard, or a sound-soft scatterer and inside a PEC core; inside a penetrable scatterer, such as a dielectric sphere, the field inside is returned.

The [`ScalarPotentialJump`](@ref dielecimped) is the difference ``\Phi_\mathrm{i} - \Phi_\mathrm{e}`` of the potential on the inner and on the outer side of the thin impedance layer.


---
## Far Fields

The far field depends on the direction of the locations only, their distance from the origin being irrelevant. The factor ``\mathrm{e}^{-\mathrm{j} k r} / r`` of the outgoing wave is omitted, that is, for the scattered electric field ``\bm e^\mathrm{sc}``
```math
\bm e^\mathrm{sc}_\infty(\hat{\bm r}) = \lim_{r \rightarrow \infty} r \, \mathrm{e}^{\mathrm{j} k r} \, \bm e^\mathrm{sc}(\bm r) \,,
```
and likewise for the acoustic pressure.

!!! note
    The total far field is not defined for a plane wave and for a spherical mode, since their incident far fields do not exist. For these excitations, `field(sp, ex, FarField(locations))` throws an error, while the scattered far field is available.


---
## Traces

The traces are evaluated on the surface of the scatterer: only the angular position of each location is taken into account, its radial coordinate being replaced by that of the surface. Hence the locations may also be given by the points of a faceted surface mesh, which do not lie exactly on the surface.

- `PressureTrace` is the Dirichlet trace ``\gamma_0 p``, the pressure on the surface.
- [`PressureNormalGradient`](@ref) is the Neumann trace ``\gamma_1 p = \hat{\bm n} \cdot \nabla p``. By default, the normal ``\hat{\bm n} = \hat{\bm r}`` of a sphere centered in the origin is employed. Normals can be passed instead, e.g., those of a surface mesh, `PressureNormalGradient(locations, normals)`; for a spheroid, `PressureNormalGradient(sp, locations)` fills in its outward normals.
- [`PressureJump`](@ref) is the jump ``[p] = p|_+ - p|_-`` across an open surface, that is, across a disc. Since the incident field is continuous across the disc, the jump of the total field equals that of the scattered one. For a closed surface an error is thrown.

!!! warning
    The outward normal of a spheroid is not parallel to the position vector, so that the default normal of [`PressureNormalGradient`](@ref) is wrong for a spheroid, see [Surface Traces](@ref).


---
## Accuracy Settings

All functions accept the keyword argument `parameter`, e.g.,
```julia
F = scatteredfield(sp, ex, quantity; parameter=Parameter(-1, 1e-8))   # nmax, relativeAccuracy
```
where the `relativeAccuracy` determines when the evaluated series are truncated, see [`Parameter`](@ref).
