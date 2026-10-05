# Spheroid and Disc

In contrast to the other pages of this section, this one is organized by scatterer rather than by excitation: the oblate and the prolate spheroid, as well as the disc, a degenerate oblate spheroid, are solved by a series in the spheroidal wave functions of their shape, which is shared by every acoustic excitation. Both the [plane wave](@ref ACpwAPI) and the [monopole](@ref ACpointAPI) are supported, at an arbitrary direction of incidence or position, and the scatterer may have an arbitrary orientation.


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
## [Solution Approach](@id ACspheroidSeries)

With the spheroidal parameter ``c = k f``, which takes the role of ``k a`` for a sphere, the incident pressure is expanded in the regular spheroidal wave functions of the shape about the center of the scatterer,
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

The oblate and prolate spheroidal wave functions are provided by [SpheroidalWaves.jl](https://github.com/Chemelli94/SpheroidalWaves.jl); for their theory see [flammerSpheroidalWaveFunctions1957](@cite). The Meixner-Schäfke normalization is employed, for which ``S_{mn}(c, \eta) \rightarrow P_n^m(\eta)`` as ``c \rightarrow 0``, so that the limit of a sphere connects directly to the series of the spherical scatterers.

#### The Coefficients of the Incident Expansion

As for the sphere, the coefficients of the incident expansion are known in closed form [flammerSpheroidalWaveFunctions1957](@cite). For a plane wave of amplitude ``a`` whose direction of incidence has the spheroidal angles ``(\eta_d, \varphi_d)`` in the frame of the scatterer,
```math
A_{mn} = 2 a \, (-\mathrm{j})^n \, \cfrac{S_{|m|n}(c, \eta_d)}{N_{|m|n}} \, \mathrm{e}^{-\mathrm{j} m \varphi_d} \,,
```
and for a monopole of amplitude ``a`` at the spheroidal position ``(\xi_s, \eta_s, \varphi_s)``, by the addition theorem of the free-space Green's function,
```math
A_{mn} = -\cfrac{\mathrm{j} k a}{2 \pi} \, \cfrac{S_{|m|n}(c, \eta_s) \, R^{(\mathrm{out})}_{|m|n}(c, \xi_s)}{N_{|m|n}} \, \mathrm{e}^{-\mathrm{j} m \varphi_s} \,,
```
valid for ``\xi < \xi_s``, which includes the surface of the scatterer. Here ``N_{mn} = \int_{-1}^{1} S_{mn}^2 \, \mathrm{d}\eta = \frac{2}{2n+1} \frac{(n+m)!}{(n-m)!}`` is the norm of the angular functions in the Meixner-Schäfke normalization, which does not depend on ``c``. As ``c \rightarrow 0`` the two reduce to ``(2n+1)(-\mathrm{j})^n`` and ``-\mathrm{j}k / (4\pi) \, (2n+1) \, h_n^{(2)}(k r_0)`` of the [plane wave](@ref ACpwAPI) and the [monopole](@ref ACpointAPI) scattered by a sphere.

!!! note
    The angular functions enter these coefficients squared, ``S(\eta_d) \, S(\eta)`` once the expansion is formed, so that their normalization and the Condon-Shortley phase cancel; only the norm ``N_{mn}`` and the normalization of the radial functions, which is the standard one, have to match.

For any other excitation, whose expansion is not known, the incident field is **projected** onto the angular functions on a surface ``\xi = \mathrm{const}``: the angular functions of equal order are orthogonal there with unit weight, and the exponentials in ``\varphi``, so that
```math
A_{mn} = \cfrac{1}{R^{(1)}_{mn}(c, \xi) \, N_{mn}} \int_{-1}^{1} g_m(\eta) \, S_{mn}(c, \eta) \, \mathrm{d}\eta \,,
\qquad
g_m(\eta) = \cfrac{1}{2\pi} \int_0^{2\pi} p_\mathrm{i} \, \mathrm{e}^{-\mathrm{j} m \varphi} \, \mathrm{d}\varphi \,.
```
The projection requires nothing but the incident field, and it does not depend on any convention of the wave functions, as the same normalization enters the integral and ``N_{mn}``, which it evaluates by the same quadrature. This makes it the fallback for new excitations, and an independent check of the closed forms above: for the plane wave and the monopole the two agree to machine precision, see [How the Implementation Is Validated](@ref ACspheroidValidation).

#### Reusing the Modal Coefficients

The spheroidal wave functions are expensive compared to the spherical ones. The coefficients are therefore computed once for a given scatterer and excitation and are reused for every evaluation point. They are also accessible directly, in order to be reused across several calls:
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
    The convenience signatures without `md` determine the coefficients internally and discard them afterwards. Whenever more than one quantity is of interest for the same scatterer and excitation, passing `md` explicitly avoids recomputing them.


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

1. **The spheroidal wave functions.** They are evaluated by [SpheroidalWaves.jl](https://github.com/Chemelli94/SpheroidalWaves.jl), whose accuracy sets a floor that this package cannot improve upon. The scattering coefficients ``b_{mn}`` are exact ratios of those functions, so their accuracy is entirely that of this layer.
2. **The coefficients ``A_{mn}`` of the incident expansion.** For the plane wave and the monopole they are closed forms of the wave functions and add no error of their own. Only for an excitation without a known expansion does the projection contribute the error of its quadrature, see [The Projection](@ref ACspheroidProjection).
3. **The truncation of the double series** at the order ``M`` and the degree ``N``.

The boundary condition, by contrast, contributes no error of its own: it is satisfied term by term, exactly, by construction of ``b_{mn}``.

#### Truncation of the Series

The degree and the order are treated differently, because they are limited by different things.

The **degree** is bounded by the size of the scatterer, as the scattering coefficients decay once a mode is cut off at the surface. The counterpart of ``k a`` is ``c \rho``, where ``\rho`` is the radius of the circumscribing sphere in units of the semifocal distance, ``\sqrt{1 + \xi_0^2}`` for an oblate and ``\xi_0`` for a prolate spheroid, so that
```math
N = \left\lceil c \rho \right\rceil + 15
```
is taken initially. The margin of 15 is a safety factor: in contrast to the spherical series, which simply adds terms until they stop contributing, the truncation has to be fixed before any field is evaluated, because the coefficients are tabulated once and reused for every evaluation point.

The **order** is not estimated but *measured*. The orders are computed in increasing ``|m|``, the size of their modes on the surface of the scatterer being measured as described below, until two consecutive ones contribute less than the relative accuracy of [`Parameter`](@ref); the last one that does is the truncation. Orders the excitation does not excite are skipped altogether: an axial excitation is rotationally symmetric and costs the order zero alone, whereas a grazing one needs ``M`` comparable to ``N``. Stopping early also keeps clear of the wave functions of high order, which exceed the range of the floating-point numbers once ``n + m`` approaches about 170.

Finally, the truncation is **verified rather than trusted**, and on the surface of the scatterer, where it matters most. The size of a mode of the scattered pressure there is its ``L^2`` norm,
```math
\left| A_{mn} \, b_{mn} \, R^{(\mathrm{out})}_{mn}(c, \xi_0) \right| \sqrt{N_{mn}} \,,
```
which bounds its contribution everywhere outside, as the outgoing radial functions decrease outward. The largest such size among the last two degrees — two, since on a disc only one parity of ``n - m`` scatters, so that the last degree alone may vanish — has to be negligible compared with the largest one overall.

!!! note
    The product ``A_{mn} b_{mn}`` alone would measure the far field only, where all outgoing radial functions are of the same size. On the surface the outgoing functions of high degree are large, and a mode with a negligible far field may contribute markedly there: for a monopole close to a disc, ``A_{mn} b_{mn}`` had decayed to ``10^{-26}`` at the cutoff while the boundary condition was violated by ``10^{-3}``.

If the omitted modes are not negligible and the degree has been determined automatically, it is **raised** by half and the coefficients are recomputed, until they are negligible, up to four times the initial degree and as long as the result keeps improving. This is what a source close to the scatterer requires, see below; for every other excitation the initial degree suffices and nothing is recomputed. A message is printed if the relative accuracy is not attained, saying whether a larger degree may help.

!!! tip
    Both truncations can be overridden. The `nmax` of [`Parameter`](@ref) fixes the degree, and the keyword arguments `M` and `N` of [`modes`](@ref SphericalScattering.modes) fix both; a given degree is used as it is, without being raised. The automatic settings have been checked against a deliberately generous truncation and agree to a relative error below ``10^{-9}``.

#### Sources Close to the Scatterer

The expansion of the field of a monopole converges the more slowly the closer the monopole is to the scatterer, so that a nearby source requires more degrees than the size of the scatterer suggests; the degree is raised automatically, see above. For a monopole above a disc of radius ``a`` and ``ka = 1.5``, for instance, the automatic degree of 17 is raised as follows:

| height of the monopole | degree | residual of the boundary condition |
|:---------------------- |:------ |:---------------------------------- |
| ``0.6 \, a`` | 59 | ``10^{-12}``, from ``10^{-5}`` |
| ``0.3 \, a`` | 68, the limit | ``10^{-10}``, from ``10^{-3}``; ``N = 89`` given explicitly attains ``2 \cdot 10^{-12}`` |
| ``0.15 \, a`` | 68, the limit | ``2 \cdot 10^{-6}``, from ``2 \cdot 10^{-2}`` |

!!! note
    The attainable accuracy is limited for a source very close to the scatterer. For the closest monopole above, the residual stalls at about ``3 \cdot 10^{-7}`` whatever the truncation, with the analytic and with the projected coefficients alike. It is not caused by cancellation in the series, whose largest term is of the size of its sum; the most likely cause is the accuracy of the wave functions for these arguments, which has not been confirmed, though. Beyond ``N \approx 130`` the wave functions exceed the range of the floating-point numbers, which [`modes`](@ref SphericalScattering.modes) reports as an error rather than returning coefficients that are not finite.

#### [The Projection](@id ACspheroidProjection)

For an excitation without a known expansion, the coefficients are obtained by projecting its field, see [`projectedCoefficients`](@ref SphericalScattering.projectedCoefficients). The surface of the projection has to lie between the scatterer and the source of the field: beyond the source the expansion in the regular wave functions does not hold, and a projection there yields coefficients that are simply wrong. By default the surface is that of the scatterer itself, or ``\xi = 0.5`` for an oblate scatterer flatter than that, and it is moved halfway between the scatterer and the source if that is closer; a given surface that reaches the source is rejected with an error. For a prolate spheroid the surface itself always serves, even for a thin one close to the line segment.

!!! note
    A disc cannot be projected on its own surface. At ``\xi = 0`` the regular radial functions of odd ``n - m`` vanish, and the projection divides by ``R^{(1)}_{mn}(c, \xi)``; the surface is therefore moved off the degenerate one.

The two integrals are discretized in ``\eta`` by Gauss-Legendre quadrature with ``2N + 16`` nodes, provided by [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), and in ``\varphi`` by ``4N + 8`` equispaced nodes, which for a periodic integrand is spectrally accurate. The norm ``N_{mn}`` is computed by the same quadrature as the integral it divides, which makes the result independent of the normalization of the angular functions.

!!! note
    The projected coefficients of high degree are accurate only as far as they matter on the surface of the projection: the error of the quadrature is amplified by ``1 / R^{(1)}_{mn}(c, \xi)``, which is large for degrees beyond ``c \xi``, whereas the field there involves the product ``A_{mn} R^{(1)}_{mn}``. Compared by the size of the modes on that surface, the projected coefficients of the plane wave and the monopole agree with the closed forms to machine precision; for a monopole close to the scatterer, whose field varies rapidly on the surface of the projection, the default quadrature limits the agreement to about ``10^{-10}``, a refined one restoring it.

#### The Outgoing Radial Function

The radial functions come in four kinds: ``R^{(1)}`` and ``R^{(2)}`` are the two real solutions, while kinds 3 and 4 are their complex combinations ``R^{(1)} \pm \mathrm{j} R^{(2)}``.

!!! warning
    The spheroidal literature calls the outgoing function ``R^{(3)}``, which corresponds to the time convention ``\mathrm{e}^{-\mathrm{j}\omega t}``. This package employs ``\mathrm{e}^{+\mathrm{j}\omega t}``, for which the outgoing solution is the analogue of ``h_n^{(2)} = j_n - \mathrm{j} y_n``, that is, **kind 4**. Picking kind 3 would yield an incoming wave.

This is not a cosmetic point. An incoming wave satisfies the boundary conditions exactly as well as an outgoing one, so no boundary-condition residual would reveal the mistake; it takes the radiation condition to do so, which is why it is tested explicitly.

#### Degenerate Points of the Coordinates

The spheroidal coordinates are not regular everywhere, and the places where they fail are precisely the geometric features of interest: the rim and the faces of a disc, and the axis, which meets every spheroid at its poles, the tips of a prolate one.

- **The rim of a disc.** The metric coefficient ``h_\xi = f |\eta|`` vanishes at ``\xi = 0``, ``\eta = 0``, so that the normal derivative ``\partial / \partial n = h_\xi^{-1} \partial / \partial \xi`` is singular there. This is the edge singularity of a disc, a property of the solution and not an artifact: the Neumann trace is genuinely unbounded at the rim, and a trace evaluated exactly at the rim is meaningless.

    The jump ``[p]`` is the better-behaved quantity, and the coordinates carry its edge behavior for free. Approaching the rim of the disc ``\xi_0 = 0`` from within at a distance ``\delta`` gives ``\eta = \sqrt{2\delta - \delta^2}``, while the modes retained by the parity split are precisely those odd in ``\eta``. Every term of the series therefore vanishes like ``\sqrt{\delta}``, which is the edge condition — no truncation of the series can violate it.
- **The axis.** At ``|\eta| = 1`` the basis vectors ``\hat{\bm e}_\eta`` and ``\hat{\bm e}_\varphi`` are undetermined, while ``h_\eta`` diverges and ``h_\varphi`` vanishes. The surface is smooth there nevertheless, and so is the field: the tangential gradient is finite and is carried by the orders ``m = \pm 1`` alone, whose angular functions vanish like ``\sqrt{1 - \eta^2}``, just as ``h_\varphi`` does. On the axis the Neumann trace is therefore evaluated as this limit, in which the dependence on ``\varphi`` cancels. This matters in practice, as a surface mesh of a spheroid usually has its vertices at the poles.
- **The faces of a disc.** They share their Cartesian coordinates, so the face cannot be recovered from a location. The one with ``\eta > 0`` is returned, and the other follows from the parity table above.

#### [How the Implementation Is Validated](@id ACspheroidValidation)

The boundary conditions are a necessary check, but on their own they are a weak one, for the reason given above: ``b_{mn}`` is *defined* by the boundary condition, so the residual reduces to the error of the incident re-expansion on the surface. The further checks close that gap, each being independent of the series under test:

| check | what it would catch |
|:----- |:------------------- |
| **Boundary conditions.** ``\gamma_0 p`` on a sound-soft and ``\gamma_1 p`` on a sound-hard scatterer, with no truncation specified. | an error in the incident coefficients or the automatic truncation |
| **Analytic against projected coefficients.** The closed forms of the plane wave and the monopole against the projection of their fields, which relies on no convention of the wave functions; compared by the size of each mode on the surface of the projection, they agree to machine precision for both shapes, the disc, and a monopole close to it. | a wrong normalization or phase of the wave functions, or a wrong closed form |
| **Sphere bridge.** As the eccentricity vanishes at fixed ``c \xi_0 = k a``, the spheroid degenerates into a sphere. The deviation from the independently validated spherical series is of the order of the eccentricity and drops by two decades when it does. | a wrong scattering coefficient, a wrong normalization, a wrong geometric parametrization |
| **Radiation condition.** ``r \, \mathrm{e}^{\mathrm{j} k r} p^\mathrm{sc}`` approaches a constant over a range of radii for the outgoing radial function, and is rejected by a wide margin for the regular one. | kind 3 in place of kind 4, that is, an incoming scattered wave |
| **Reference truncation.** The automatic settings against a deliberately generous one. | an insufficient automatic truncation |
| **The axis.** The Neumann trace on the axis against the average over two antipodal points next to it, which approaches the axis quadratically in their distance. | a wrong limit of the tangential gradient |

Underneath, the wave-function layer is checked on its own terms: the Wronskian ``R^{(1)} R^{(2)\prime} - R^{(2)} R^{(1)\prime} = 1 / (c (\xi^2 \pm 1))``, the upper sign for the oblate and the lower for the prolate shape, the limit ``S_{mn}(c, \eta) \rightarrow P_n^m(\eta)`` as ``c \rightarrow 0``, the asymptotics of the outgoing function, and the parity split of the scattering coefficients on a disc.

All of these checks are carried out for both shapes, the prolate spheroid being tested in addition close to its excluded degeneration, for ``a / b = 0.05``, and at its tips, where the curvature is largest.


---
## Limitations

!!! note
    - A sphere is not a limiting case that can be evaluated: the coordinates degenerate, hence ``a \neq b`` is enforced. Use [`HardSphere`](@ref) or [`SoftSphere`](@ref).
    - The prolate spheroid degenerating into a line segment, ``a = 0``, is not supported; a thin one close to it is, and is as accurate as any other.
    - The jump of the **normal derivative** across a disc, the natural unknown of the sound-soft case, is not implemented.
    - A monopole has to lie outside the scatterer, as the expansion of its field assumes. This is checked, an error being thrown otherwise.
