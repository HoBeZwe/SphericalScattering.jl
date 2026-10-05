# Spheroid and Disc

The oblate and the prolate spheroid, as well as the disc, a degenerate oblate spheroid, are solved by a series in the spheroidal wave functions of their shape, which is shared by every acoustic excitation. Both the [plane wave](@ref ACpwAPI) and the [monopole](@ref ACpointAPI) are supported, at an arbitrary direction of incidence or position, and the scatterer may have an arbitrary orientation. How the solution is obtained, and how accurate it is, is described under [Spheroidal Solution](@ref ACspheroidSolution).


---
## Geometry

The scatterer is a surface ``\xi = \xi_0`` of spheroidal coordinates ``(\xi, \eta, \varphi)``. It is specified by its equatorial radius ``a`` and its polar radius ``b``, the latter measured along the axis of revolution, and the shape follows from them:

| | radii | semifocal distance ``f`` | surface | coordinates |
|:-- |:-- |:-- |:-- |:-- |
| [`OblateSpheroid`](@ref) | ``a > b`` | ``\sqrt{a^2 - b^2}`` | ``\xi_0 = b / f \ge 0`` | [oblate](@ref oblateCoords) |
| [`ProlateSpheroid`](@ref) | ``b > a`` | ``\sqrt{b^2 - a^2}`` | ``\xi_0 = b / f > 1`` | [prolate](@ref prolateCoords) |

The special case ``b = 0`` of the oblate spheroid, that is, ``\xi_0 = 0``, is the disc of radius ``a = f``, whose two faces are distinguished by the sign of ``\eta``. The corresponding degeneration of the prolate spheroid, ``a = 0`` or ``\xi_0 = 1``, is the segment of the axis between the foci, which is not supported.

!!! note
    ``a \neq b`` is required: the spheroidal coordinates degenerate for a sphere, where ``f \rightarrow 0`` and ``\xi_0 \rightarrow \infty``. A sphere is a [`HardSphere`](@ref) or a [`SoftSphere`](@ref), for which the considerably cheaper spherical series is evaluated.

The boundary condition is carried as a type parameter, as it is orthogonal to the geometry:

| type | condition on the surface |
|:---- |:------------------------ |
| [`SoundHard`](@ref) | ``\gamma_1 p = \hat{\bm n} \cdot \nabla p = 0`` (the normal velocity vanishes) |
| [`SoundSoft`](@ref) | ``\gamma_0 p = 0`` (pressure release) |

#### [API](@id ACspheroidAPI)

```@docs
Spheroid
Disc
```

The axis of revolution ``\hat{\bm a}`` is normalized during the initialization and defaults to ``\hat{\bm e}_z``. Since the pressure is a scalar, an arbitrarily oriented spheroid costs no more than transforming the evaluation points into the frame of the scatterer; nothing has to be rotated back.

The dimensions can be queried, which is of interest mainly because a `Spheroid` stores ``f`` and ``\xi_0`` rather than the radii it was constructed from:
```@docs
equatorialRadius
polarRadius
isdisc
```


---
## Quantities

The general API is employed, with or without the modal coefficients:
```julia
p  = scatteredfield(sp, ex, Pressure(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))

γ₀ = scatteredfield(sp, ex, PressureTrace(point_cart))

γ₁ = scatteredfield(sp, ex, PressureNormalGradient(sp, point_cart))

Δp = scatteredfield(sp, ex, PressureJump(point_cart))     # a disc only
```
and likewise `field` for the total quantities, where `sp` is a [`Spheroid`](@ref) — possibly the degenerate one returned by [`Disc`](@ref). As for a sphere, the total far field exists for a monopole but not for a plane wave, and the pressure is set to zero inside the scatterer, where the oblate spheroidal coordinates answer what the radius of a sphere cannot.

The far field omits the factor ``\mathrm{e}^{-\mathrm{j} k r} / r``, as for the other scatterers, see [Quantities](@ref quantitiesConcept); how the quantities follow from the modes is described under [The Quantities in Terms of the Modes](@ref ACspheroidModeQuantities).

#### Surface Traces

The traces are evaluated on the surface: only the direction of each location is taken into account, its radial coordinate being replaced by ``\xi_0``. Hence the locations may also be given by the points of a faceted surface mesh, which do not lie exactly on the spheroid.

!!! warning
    The outward normal of a spheroid is **not** parallel to the position vector — for a disc it is ``\pm \hat{\bm e}_z`` everywhere. The default of [`PressureNormalGradient`](@ref), the radial direction, is therefore wrong for a spheroid. Pass the scatterer, which fills in the correct normals:
    ```julia
    γ₁ = scatteredfield(sp, ex, PressureNormalGradient(sp, point_cart))
    ```
```@docs
outwardNormal
outwardNormals
```

#### [The Jump Across a Disc](@id ACdiscJump)

On an open surface a one-sided trace is not determined by the geometry alone, which is why the natural unknown of a boundary element formulation is the jump ``[p] = p|_+ - p|_-``, see [`PressureJump`](@ref). The two faces of the degenerate surface ``\xi = 0`` are ``\eta > 0`` and ``\eta < 0``. The jump of the total field equals that of the scattered field, the incident field being continuous across the disc.

!!! note
    At the degenerate surface the radial functions split by parity, and the two boundary conditions select opposite halves of the spectrum: the Neumann problem is carried by the modes with odd ``n - m``, the Dirichlet problem by those with even ``n - m``. Consequently ``[p]`` is carried entirely by the sound-hard disc and vanishes identically for the sound-soft one, whose unknown is the jump of the normal derivative instead — that quantity is not implemented yet.

The same parity governs the one-sided traces, for which the face with ``\eta > 0`` is returned:

| disc | ``\gamma_0`` | ``\gamma_1`` |
|:---- |:------------ |:------------ |
| sound-soft | same on both faces | opposite |
| sound-hard | opposite | vanishes |


---
## Reusing the Modal Coefficients

The spheroidal wave functions are expensive compared to the spherical ones. The coefficients are therefore computed once for a given scatterer and excitation and are reused for every evaluation point. They are also accessible directly, in order to be reused across several calls:
```julia
md = SphericalScattering.modes(sp, ex)                 # the expensive part

p  = scatteredfield(sp, ex, md, Pressure(points))      # cheap from here on
FF = scatteredfield(sp, ex, md, FarField(directions))
```
```@docs
SphericalScattering.modes
```

What is stored is the spheroidal parameter ``c``, the two truncations, and the two coefficient tables ``A_{mn}`` and ``b_{mn}``, indexed as `[m + M + 1, n + 1]`, see [Solution Approach](@ref ACspheroidSeries) for their meaning, and see [`SpheroidalModes`](@ref SphericalScattering.SpheroidalModes).

!!! tip
    The convenience signatures without `md` determine the coefficients internally and discard them afterwards. Whenever more than one quantity is of interest for the same scatterer and excitation, passing `md` explicitly avoids recomputing them.


---
## Limitations

!!! note
    - A sphere is not a limiting case that can be evaluated: the coordinates degenerate, hence ``a \neq b`` is enforced. Use [`HardSphere`](@ref) or [`SoftSphere`](@ref).
    - The prolate spheroid degenerating into a line segment, ``a = 0``, is not supported; a thin one close to it is, and is as accurate as any other.
    - The jump of the **normal derivative** across a disc, the natural unknown of the sound-soft case, is not implemented.
    - A monopole has to lie outside the scatterer, as the expansion of its field assumes. This is checked, an error being thrown otherwise.
