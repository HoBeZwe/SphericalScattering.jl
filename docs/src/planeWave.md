
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

The plane wave with amplitude ``a``, wave vector ``\bm k``, and polarization ``\hat{\bm p}`` is assumed to have the field
```math
\bm e_\mathrm{PW}(\bm r) = a \hat{\bm p}  \, \mathrm{e}^{-\mathrm{j} \bm k \cdot \bm r}  \,.
```


---
## [API](@id pwAPI)

The API provides the following constructor with default values:
```@docs
planeWave
```

---
## Incident Field

The electric field of the plane wave is as given above. The magnetic field is given by
```math
\bm h_\mathrm{PW}(\bm r) = \cfrac{a}{Z_\mathrm{F}} (\hat{\bm k} \times \hat{\bm p})  \mathrm{e}^{-\mathrm{j} \bm k \cdot \bm r}
```
with ``Z_\mathrm{F} = \sqrt{\mu / \varepsilon}``

#### API

The general API is employed:
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
```

---
## Scattered Field

The scattered field computation follows [[1, pp. 347ff]](@ref refs). 

!!! warning
    So far the plane wave is assumed to travel in positive ``z``-axis direction and to have a polarization along the ``x``-axis! This is planned to be generalized.

#### API

The general API is employed:
```julia
E  = scatteredfield(sp, ex, ElectricField(point_cart))

H  = scatteredfield(sp, ex, MagneticField(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))
```

---
## Total Field

!!! warning
    Not fully implemented yet.

#### API

The general API is employed:
```julia
E  = field(sp, ex, ElectricField(point_cart))

H  = field(sp, ex, MagneticField(point_cart))

FF = field(sp, ex, FarField(point_cart))
```