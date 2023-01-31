
# Plane Wave

```@raw html
<figure>
  <img
    src="../assets/PECsphere.svg"
    alt="Setup"
    width="300" />

  <figcaption>
    Setup for the excitation with a plane wave. TODO: more suitable image.
  </figcaption>
</figure>
<br/>
```

---
## Definition

A plane wave with amplitude ``a``, wave vector ``\bm k = k \hat{\bm k}``, and polarization ``\hat{\bm p}`` (vectors with a hat denote unit vectors) is defined by the field
```math
\bm e_\mathrm{PW}(\bm r) = a \hat{\bm p}  \, \mathrm{e}^{-\mathrm{j} \bm k \cdot \bm r}  \,,
```
where the ploarization and wave vector are orthogonal, that is,
```math
\bm k \cdot \hat{\bm p} = 0
```
holds.

---
## [API](@id pwAPI)

The API provides the following constructor with default values:
```@docs
planeWave
```

!!! note
    The provided `direction` and the `polarization` vectors have to be orthogonal. This is checked during initialization.

!!! tip
    The `direction` and the `polarization` vectors are each automatically normalized to unit vectors during the initialization.

---
## Incident Field

The electric field of the plane wave is as given above. The magnetic field is given by
```math
\bm h_\mathrm{PW}(\bm r) = \cfrac{a}{Z_\mathrm{F}} (\hat{\bm k} \times \hat{\bm p})  \mathrm{e}^{-\mathrm{j} \bm k \cdot \bm r}
```
with ``Z_\mathrm{F} = \sqrt{\mu / \varepsilon}``.

#### API

The general API is employed:
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))
```

!!! note
    The far-field of a plane wave is not defined.

---
## Scattered Field

The scattered field computation follows [[1, pp. 347ff]](@ref refs). 

!!! note
    Internal details of the computations: Following [[1, pp. 347ff]](@ref refs) the plane wave is initially assumed to travel in positive ``z``-axis direction and to have a polarization along the positive ``x``-axis. Arbitrary directions and orientations (forming a valid pair) are obtained via [rotations](@ref rotationDetails). 

#### API

The general API is employed:
```julia
E  = scatteredfield(sp, ex, ElectricField(point_cart))

H  = scatteredfield(sp, ex, MagneticField(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))
```

---
## Total Field

#### API

The general API is employed:
```julia
E  = field(sp, ex, ElectricField(point_cart))

H  = field(sp, ex, MagneticField(point_cart))
```

!!! note
    The total far-field is not defined (since the incident far-field is not defined).