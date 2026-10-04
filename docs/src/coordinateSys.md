
# Coordinate System

!!! info
    At the user interface a Cartesian basis is used for all coordinates and field components. Internally, spherical, oblate spheroidal, and prolate spheroidal coordinates are used as well.


## Spherical Coordinates

The employed coordinate system uses the following convention for the spherical coordinates everywhere in the code with ``\vartheta \in [0, \pi]`` and ``\varphi \in [-\pi, \pi]``.

```@raw html
<div align="center">
<img src="../assets/CoordinateSystem.svg" width="350"/>
</div>
<br/>
```

```@raw html
<!---
<center>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/SphereEx_white.svg" height="3000" align="middle">
  <source media="(prefers-color-scheme: light)" srcset="assets/SphereEx.svg" height="3000">
  <img alt="" src="" height="3000">
</picture>
</center>
--->
```


## [Oblate Spheroidal Coordinates](@id oblateCoords)

The oblate spheroidal coordinates ``(\xi, \eta, \varphi)`` with ``\xi \ge 0``, ``\eta \in [-1, 1]`` and ``\varphi \in [-\pi, \pi]`` are defined by
```math
x = f \sqrt{(1 + \xi^2)(1 - \eta^2)} \cos \varphi \,, \qquad
y = f \sqrt{(1 + \xi^2)(1 - \eta^2)} \sin \varphi \,, \qquad
z = f \xi \eta
```
where ``f`` denotes the semifocal distance, so that the distance from the origin is ``r = f \sqrt{\xi^2 + 1 - \eta^2}``.

The surface ``\xi = \xi_0`` is an oblate spheroid with equatorial radius ``f \sqrt{1 + \xi_0^2}`` and polar radius ``f \xi_0``; the surfaces ``\eta = \mathrm{const}`` are the confocal one-sheeted hyperboloids. The coordinate ``\eta`` thus plays the role of ``\cos \vartheta`` and reduces to it in the limit of large ``\xi``, whereas ``\varphi`` is the azimuth of the spherical coordinates.

Two degenerate cases are of interest. For ``\xi_0 \rightarrow \infty`` the spheroid approaches a sphere, with ``f \rightarrow 0`` such that ``f \xi_0`` stays finite, and the coordinates degenerate. For ``\xi_0 = 0`` the surface collapses onto the disc of radius ``f`` in the plane ``z = 0``, its two faces being distinguished by the sign of ``\eta``: inside the disc, the sign of ``\eta`` cannot be recovered from the Cartesian coordinates of a point.

The metric coefficients read
```math
h_\xi = f \sqrt{\cfrac{\xi^2 + \eta^2}{1 + \xi^2}} \,, \qquad
h_\eta = f \sqrt{\cfrac{\xi^2 + \eta^2}{1 - \eta^2}} \,, \qquad
h_\varphi = f \sqrt{(1 + \xi^2)(1 - \eta^2)} \,.
```
Since ``\hat{\bm e}_\xi = h_\xi^{-1} \partial \bm r / \partial \xi`` is the outward normal of the surface ``\xi = \mathrm{const}``, the normal derivative is ``\partial / \partial n = h_\xi^{-1} \partial / \partial \xi``.

!!! note
    The basis degenerates at the rim of a disc, where ``\xi = 0`` and ``\eta = 0``, and on the axis, where ``|\eta| = 1``. These are the edge and the poles of the scatterer; see the [accuracy](@ref ACspheroidAccuracy) of the spheroidal solution for what follows from it.

!!! warning
    The outward normal ``\hat{\bm e}_\xi`` is not parallel to the position vector unless the scatterer is a sphere. For a disc it is ``\pm \hat{\bm e}_z`` everywhere.

These coordinates are employed for the [acoustic spheroid and disc](@ref ACspheroidAPI). The conversions are available as `SphericalScattering.obl2cart`, `SphericalScattering.cart2obl`, `SphericalScattering.oblateMetric` and `SphericalScattering.oblateBasis`.


## Prolate Spheroidal Coordinates

Text