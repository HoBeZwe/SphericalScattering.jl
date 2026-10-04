# Spheroid and Disc

In contrast to the other pages of this section, this one is organized by scatterer rather than by excitation: the oblate spheroid and its degenerate case, the disc, are solved by a series in the oblate spheroidal wave functions, which is shared by every acoustic excitation. Both the [plane wave](@ref ACpwAPI) and the [monopole](@ref ACpointAPI) are supported, at an arbitrary direction of incidence or position, and the scatterer may have an arbitrary orientation.


---
## Geometry

The scatterer is a surface ``\xi = \xi_0`` of the oblate spheroidal coordinates ``(\xi, \eta, \varphi)``, see the [coordinate system](@ref oblateCoords). It is specified by its equatorial radius ``a`` and its polar radius ``b < a``, from which
```math
f = \sqrt{a^2 - b^2} \qquad \text{and} \qquad \xi_0 = \cfrac{b}{f}
```
follow, the semifocal distance ``f`` and the radial coordinate of the surface. The special case ``b = 0``, that is, ``\xi_0 = 0``, is the disc of radius ``a = f``, whose two faces are distinguished by the sign of ``\eta``.

!!! note
    Strictly ``a > b`` is required: the oblate spheroidal coordinates degenerate for a sphere, where ``f \rightarrow 0`` and ``\xi_0 \rightarrow \infty``. A sphere is a [`HardSphere`](@ref) or a [`SoftSphere`](@ref), for which the considerably cheaper spherical series is evaluated.

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
## [Solution Approach](@id ACspheroidSeries)

With the spheroidal parameter ``c = k f``, which takes the role of ``k a`` for a sphere, the incident pressure is expanded in the regular oblate spheroidal wave functions about the center of the scatterer,
```math
p_\mathrm{i}(\bm r) = \sum_{m=-M}^{M} \sum_{n=|m|}^{N} A_{mn} R^{(1)}_{mn}(c, \xi) \, S_{mn}(c, \eta) \, \mathrm{e}^{\mathrm{j} m \varphi} \,,
```
where ``R^{(1)}_{mn}`` and ``S_{mn}`` denote the radial and the angular function, the counterparts of ``j_n`` and ``P_n^m``. Since ``\xi = \xi_0`` is a coordinate surface and the angular functions form a complete orthogonal set on it, the boundary condition decouples the modes: it is satisfied term by term, so that the scattered pressure follows from the incident expansion by multiplying each coefficient with the scattering coefficient of its mode and by replacing the regular radial function with the outgoing one,
```math
p^\mathrm{sc}(\bm r) = \sum_{mn} A_{mn} \, b_{mn} \, R^{(\mathrm{out})}_{mn}(c, \xi) \, S_{mn}(c, \eta) \, \mathrm{e}^{\mathrm{j} m \varphi} \,.
```
Applying the two boundary conditions yields
```math
b_{mn} = -\cfrac{R^{(1)\prime}_{mn}(c, \xi_0)}{R^{(\mathrm{out})\prime}_{mn}(c, \xi_0)}
\qquad \text{and} \qquad
b_{mn} = -\cfrac{R^{(1)}_{mn}(c, \xi_0)}{R^{(\mathrm{out})}_{mn}(c, \xi_0)}
```
for a sound-hard and a sound-soft spheroid, respectively, the prime denoting the derivative with respect to ``\xi`` [bowmanElectromagneticAcousticScattering1970](@cite). As for a sphere, the scattering coefficients do not depend on the excitation: only the coefficients ``A_{mn}`` of the incident expansion do.

The oblate spheroidal wave functions are provided by [SpheroidalWaves.jl](https://github.com/Chemelli94/SpheroidalWaves.jl); for their theory see [flammerSpheroidalWaveFunctions1957](@cite). The Meixner-Schäfke normalization is employed, for which ``S_{mn}(c, \eta) \rightarrow P_n^m(\eta)`` as ``c \rightarrow 0``, so that the limit of a sphere connects directly to the series of the spherical scatterers.

#### The Coefficients of the Incident Expansion

In contrast to the spherical case, no closed-form expansion coefficients are employed. Instead the incident field, which is known in closed form for every excitation, is **projected** onto the angular functions on a surface ``\xi = \mathrm{const}``. The angular functions of equal order are orthogonal there with unit weight, and the exponentials are orthogonal in ``\varphi``, hence
```math
A_{mn} = \cfrac{1}{R^{(1)}_{mn}(c, \xi) \, N_{mn}} \int_{-1}^{1} g_m(\eta) \, S_{mn}(c, \eta) \, \mathrm{d}\eta \,,
\qquad
g_m(\eta) = \cfrac{1}{2\pi} \int_0^{2\pi} p_\mathrm{i} \, \mathrm{e}^{-\mathrm{j} m \varphi} \, \mathrm{d}\varphi
```
with the norm ``N_{mn} = \int_{-1}^{1} S_{mn}^2 \, \mathrm{d}\eta``.

!!! note
    The projection is deliberately agnostic: it makes no assumption about the excitation beyond its field being regular in the region of the projection surface, and it is independent of the normalization of the angular functions, since the same normalization factor enters the integral and ``N_{mn}``. A new acoustic excitation therefore requires nothing but its incident field in order to be scattered by a spheroid.

#### Reusing the Modal Coefficients

The oblate spheroidal wave functions are expensive compared to the spherical ones. The coefficients are therefore computed once for a given scatterer and excitation and are reused for every evaluation point. They are also accessible directly, in order to be reused across several calls:
```julia
md = SphericalScattering.modes(sp, ex)                 # the expensive part

p  = scatteredfield(sp, ex, md, Pressure(points))      # cheap from here on
FF = scatteredfield(sp, ex, md, FarField(directions))
```
```@docs
SphericalScattering.modes
```

What is stored is the spheroidal parameter ``c``, the two truncations, and the two coefficient tables ``A_{mn}`` and ``b_{mn}``, indexed as `[m + M + 1, n + 1]`; see [`SpheroidalModes`](@ref SphericalScattering.SpheroidalModes).

!!! tip
    The convenience signatures without `md` determine the coefficients internally and discard them afterwards. Whenever more than one quantity is of interest for the same scatterer and excitation, passing `md` explicitly avoids recomputing by far the most expensive part.


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

The far field omits the factor ``\mathrm{e}^{-\mathrm{j} k r} / r``, as for the electromagnetic excitations. Since ``R^{(\mathrm{out})}_{mn}(c, \xi) \rightarrow \mathrm{j}^{n+1} \mathrm{e}^{-\mathrm{j} c \xi} / (c \xi)`` for ``\xi \rightarrow \infty``, just as ``h_n^{(2)}`` does, and since ``r = f \sqrt{\xi^2 + 1 - \eta^2} \rightarrow f \xi``, it reads
```math
p^\mathrm{sc}_\infty(\hat{\bm r}) = \lim_{r \rightarrow \infty} r \, \mathrm{e}^{\mathrm{j} k r} p^\mathrm{sc}(\bm r)
    = \cfrac{1}{k} \sum_{mn} A_{mn} \, b_{mn} \, \mathrm{j}^{n+1} S_{mn}(c, \eta) \, \mathrm{e}^{\mathrm{j} m \varphi} \,,
```
that is, no radial function has to be evaluated at all. In the limit ``\eta`` becomes the cosine of the angle between the axis of the scatterer and the direction of observation, which is why it is taken from that direction rather than from the oblate spheroidal coordinates of the location.

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

For an arbitrary normal the gradient does not reduce to two terms, as it does for a sphere: the scattered field of a spheroid depends on ``\varphi`` as well, so that all three components contribute,
```math
\hat{\bm n} \cdot \nabla p^\mathrm{sc} = \cfrac{n_\xi}{h_\xi} \cfrac{\partial p^\mathrm{sc}}{\partial \xi}
    + \cfrac{n_\eta}{h_\eta} \cfrac{\partial p^\mathrm{sc}}{\partial \eta}
    + \cfrac{n_\varphi}{h_\varphi} \cfrac{\partial p^\mathrm{sc}}{\partial \varphi} \,,
```
where the derivative with respect to ``\varphi`` amounts to a factor ``\mathrm{j}m``. All three are included, so that providing normals is supported here as well; for the outward normal ``\hat{\bm n} = \hat{\bm e}_\xi`` the tangential contributions drop out.

#### The Jump Across a Disc

On an open surface a one-sided trace is not determined by the geometry alone, which is why the natural unknown of a boundary element formulation is the jump. The two faces of the degenerate surface ``\xi = 0`` are ``\eta > 0`` and ``\eta < 0``, and the angular functions obey ``S_{mn}(c, -\eta) = (-1)^{n-m} S_{mn}(c, \eta)``, so that
```math
[p] = p|_+ - p|_- = 2 \sum_{n - m \ \mathrm{odd}} A_{mn} \, b_{mn} \, R^{(\mathrm{out})}_{mn}(c, 0) \, S_{mn}(c, \eta) \, \mathrm{e}^{\mathrm{j} m \varphi} \,,
```
the modes of even ``n - m`` cancelling. Weighting the terms by their parity avoids having to evaluate the two faces separately, which the Cartesian coordinates of a location cannot distinguish. The jump of the total field equals this one, the incident field being continuous across the disc, see [`PressureJump`](@ref).

!!! note
    At the degenerate surface the radial functions split by parity, and the two boundary conditions select opposite halves of the spectrum: the Neumann problem is carried by the modes with odd ``n - m``, the Dirichlet problem by those with even ``n - m``. Consequently ``[p]`` is carried entirely by the sound-hard disc and vanishes identically for the sound-soft one, whose unknown is the jump of the normal derivative instead — that quantity is not implemented yet.

The same parity governs the one-sided traces, for which the face with ``\eta > 0`` is returned:

| disc | ``\gamma_0`` | ``\gamma_1`` |
|:---- |:------------ |:------------ |
| sound-soft | same on both faces | opposite |
| sound-hard | opposite | vanishes |


---
## [Accuracy](@id ACspheroidAccuracy)

A spheroidal solution is a good deal more delicate than a spherical one, and most of what follows has no counterpart on the other acoustic pages. Where the error comes from, and what the package does about it, is therefore spelled out here.

#### Sources of Error

There are three, and they are independent of one another:

1. **The oblate spheroidal wave functions.** They are evaluated by [SpheroidalWaves.jl](https://github.com/Chemelli94/SpheroidalWaves.jl), whose accuracy sets a floor that this package cannot improve upon. The scattering coefficients ``b_{mn}`` are exact ratios of those functions, so their accuracy is entirely that of this layer.
2. **The projection that yields the coefficients ``A_{mn}``.** This contributes a quadrature error and, more importantly, is subject to the conditioning of the projection surface discussed below.
3. **The truncation of the double series** at the order ``M`` and the degree ``N``.

The boundary condition, by contrast, contributes no error of its own: it is satisfied term by term, exactly, by construction of ``b_{mn}``.

#### Truncation of the Series

The degree and the order are treated differently, because they are limited by different things.

The **degree** is bounded by the size of the scatterer, as the scattering coefficients decay once a mode is cut off at the surface. The counterpart of ``k a`` is ``c \sqrt{1 + \xi_0^2}``, the radicand being the equatorial radius in units of the semifocal distance, so that
```math
N = \left\lceil c \sqrt{1 + \xi_0^2} \right\rceil + 15
```
is taken. The margin of 15 is a safety factor: in contrast to the spherical series, which simply adds terms until they stop contributing, the truncation has to be fixed in advance here, because the coefficients are tabulated once and reused for every evaluation point.

The **order** is not estimated but *measured*. The azimuthal spectrum ``g_m`` of the incident field on the projection surface is evaluated for all ``|m| \le N``, and the orders are retained as long as they contribute more than the relative accuracy of [`Parameter`](@ref). This matters in practice, since every order costs a pair of calls to the spheroidal wave functions: an axial excitation is rotationally symmetric and needs ``M = 0``, whereas a grazing one needs ``M`` comparable to ``N``.

Finally, the truncation is **verified rather than trusted**. The contribution of a mode to the scattered field is carried by the product ``A_{mn} b_{mn}``, which has to have decayed at the cutoff. After assembling the coefficients, the largest contribution at ``n = N`` is compared with the largest one overall, and a message is printed should the tail not have decayed below the relative accuracy. Seeing that message is the signal to raise the truncation by hand.

!!! tip
    Both truncations can be overridden. The `nmax` of [`Parameter`](@ref) fixes the degree, and the keyword arguments `M` and `N` of [`modes`](@ref SphericalScattering.modes) fix both. The automatic settings have been checked against a deliberately generous truncation and agree to a relative error below ``10^{-9}``.

#### The Projection Surface

The surface on which the incident field is projected has to lie within the region in which that field is regular, that is, closer to the scatterer than a monopole. By default it is the surface of the scatterer itself, or ``\xi = 0.5`` for a scatterer flatter than that.

!!! note
    A disc cannot be projected on its own surface. At ``\xi = 0`` the regular radial functions of odd ``n - m`` vanish, and the projection divides by ``R^{(1)}_{mn}(c, \xi)``; the surface is therefore moved off the degenerate one. This is also why the default is a maximum rather than simply ``\xi_0``.

!!! warning
    The conditioning of the projection concerns the **incident re-expansion alone**, not the scattered field. Reconstructing ``p_\mathrm{i}`` from the coefficients ``A_{mn}`` at a radial coordinate much larger than the projection surface amplifies the error, because the regular radial functions grow steeply with ``\xi`` at a fixed degree. For that purpose — and only for it — the keyword `ξ` should be chosen comparable to the largest radial coordinate of interest, and `N` large enough that ``N \gtrsim c \xi``. The scattered field is not affected: there the outgoing radial functions and the scattering coefficients both decay with the degree, so that a truncation sufficient near the scatterer is sufficient everywhere.

#### Quadrature

The two integrals of the projection are discretized to be essentially exact for the retained modes:

- in ``\eta`` by Gauss-Legendre quadrature with ``2N + 16`` nodes, provided by [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl). The weight is unity, the angular functions being orthogonal on ``[-1, 1]`` as they stand.
- in ``\varphi`` by ``4N + 8`` equispaced nodes, which for a periodic integrand is spectrally accurate. The number is deliberately independent of ``M``, so that ``M`` can be measured from the resulting spectrum rather than having to be assumed before it.

The norm ``N_{mn}`` is computed by the same ``\eta`` quadrature as the integral it divides, which is what makes the result independent of the normalization of the angular functions.

#### The Outgoing Radial Function

The radial functions come in four kinds: ``R^{(1)}`` and ``R^{(2)}`` are the two real solutions, while kinds 3 and 4 are their complex combinations ``R^{(1)} \pm \mathrm{j} R^{(2)}``.

!!! warning
    The spheroidal literature calls the outgoing function ``R^{(3)}``, which corresponds to the time convention ``\mathrm{e}^{-\mathrm{j}\omega t}``. This package employs ``\mathrm{e}^{+\mathrm{j}\omega t}``, for which the outgoing solution is the analogue of ``h_n^{(2)} = j_n - \mathrm{j} y_n``, that is, **kind 4**. Picking kind 3 would yield an incoming wave.

This is not a cosmetic point. An incoming wave satisfies the boundary conditions exactly as well as an outgoing one, so no boundary-condition residual would reveal the mistake; it takes the radiation condition to do so, which is why it is tested explicitly.

#### Degenerate Points of the Coordinates

The oblate spheroidal coordinates are not regular everywhere, and the places where they fail are precisely the geometric features of interest.

- **The rim of a disc.** The metric coefficient ``h_\xi = f |\eta|`` vanishes at ``\xi = 0``, ``\eta = 0``, so that the normal derivative ``\partial / \partial n = h_\xi^{-1} \partial / \partial \xi`` is singular there. This is the edge singularity of a disc, a property of the solution and not an artifact: the Neumann trace is genuinely unbounded at the rim, and a trace evaluated exactly at the rim is meaningless.

    The jump ``[p]`` is the better-behaved quantity, and the coordinates carry its edge behavior for free. Approaching the rim of the disc ``\xi_0 = 0`` from within at a distance ``\delta`` gives ``\eta = \sqrt{2\delta - \delta^2}``, while the modes retained by the parity split are precisely those odd in ``\eta``. Every term of the series therefore vanishes like ``\sqrt{\delta}``, which is the edge condition — no truncation of the series can violate it.
- **The axis.** At ``|\eta| = 1`` the azimuthal direction ``\hat{\bm e}_\varphi`` is undetermined. The azimuthal derivative is set to zero there, which is exact rather than a patch: the angular functions of non-vanishing order vanish on the axis.
- **The faces of a disc.** They share their Cartesian coordinates, so the face cannot be recovered from a location. The one with ``\eta > 0`` is returned, and the other follows from the parity table above.

#### How the Implementation Is Validated

The boundary conditions are a necessary check, but on their own they are a weak one, for the reason given above: ``b_{mn}`` is *defined* by the boundary condition, so the residual reduces to the error of the incident re-expansion on the surface. Three further checks close that gap, each independent of the series under test, so that four altogether are in place:

| check | what it would catch |
|:----- |:------------------- |
| **Boundary conditions.** ``\gamma_0 p`` on a sound-soft and ``\gamma_1 p`` on a sound-hard scatterer, with no truncation specified. | an error in the projection, the quadrature or the automatic truncation |
| **Sphere bridge.** As the eccentricity vanishes at fixed ``c \xi_0 = k a``, the spheroid degenerates into a sphere. The deviation from the independently validated spherical series is of the order of the eccentricity and drops by two decades when it does. | a wrong scattering coefficient, a wrong normalization, a wrong geometric parametrization |
| **Radiation condition.** ``r \, \mathrm{e}^{\mathrm{j} k r} p^\mathrm{sc}`` approaches a constant over a range of radii for the outgoing radial function, and is rejected by a wide margin for the regular one. | kind 3 in place of kind 4, that is, an incoming scattered wave |
| **Reference truncation.** The automatic settings against a deliberately generous one. | an insufficient automatic truncation |

Underneath, the wave-function layer is checked on its own terms: the Wronskian ``R^{(1)} R^{(2)\prime} - R^{(2)} R^{(1)\prime} = 1 / (c (\xi^2 + 1))``, the limit ``S_{mn}(c, \eta) \rightarrow P_n^m(\eta)`` as ``c \rightarrow 0``, the asymptotics of the outgoing function, and the parity split of the scattering coefficients on a disc.


---
## Limitations

!!! note
    - Only the **oblate** spheroid is implemented. The prolate spheroid, for which `SpheroidalWaves.jl` is equally capable, is not yet wired up.
    - A sphere is not a limiting case that can be evaluated: the coordinates degenerate, hence ``a > b`` is enforced. Use [`HardSphere`](@ref) or [`SoftSphere`](@ref).
    - The jump of the **normal derivative** across a disc, the natural unknown of the sound-soft case, is not implemented.
    - A monopole has to lie outside the scatterer, as the expansion of its field assumes. This is checked, an error being thrown otherwise.
