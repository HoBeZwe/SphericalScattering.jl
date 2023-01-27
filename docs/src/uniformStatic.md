
# Uniform Static Field

```@raw html
<figure>
  <img
    src="../assets/PECsphere.svg"
    alt="Setup"
    width="300" />

  <figcaption>
    Setup for the excitation with a uniform static field. TODO: more suitable image.
  </figcaption>
</figure>
<br/>
```

---
## Definition

The static electric field with amplitude ``a`` and direction ``\hat{\bm p}`` is given by
```math
\bm e(\bm r) = a \, \hat{ \bm p} \,.
```


---
## [API](@id uniformAPI)

The API provides the following constructor with default values:
```@docs
UniformField
```

!!! warning
    The direction is so far not automatically normalized!


---
## Incident Field

The electric field of the plane wave is as given above. 

#### API

The general API is employed:
```julia
Φ = field(ex, ScalarPotential(point_cart))

E = field(ex, ElectricField(point_cart))
```

---
## Scattered Field

The scattered field computation follows [[3]](@ref refs). 

!!! warning
    So far the static electric field is assumed to point in ``x``-direction. This is planned to be generalized.

#### API

The following API is employed:
```julia
Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))

E = scatteredfield(sp, ex, ElectricField(point_cart))
```

---
## Total Field

!!! warning
    Not fully implemented yet.

#### API

The following API is employed:
```julia
Φ = field(sp, ex, ScalarPotential(point_cart))

E = field(sp, ex, ElectricField(point_cart))
```