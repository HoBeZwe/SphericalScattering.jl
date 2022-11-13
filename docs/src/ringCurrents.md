
# Ring Currents

```@raw html
<figure>
  <img
    src="assets/Fig_SphereRC.svg"
    alt="Setup"
    width="300" />

  <figcaption>
    Setup for the excitation with an electric or a magnetic ring current.
  </figcaption>
</figure>
<br/>
```

---
## Definition

The ring currents are defined by a circular loop of radius ``a`` carrying a uniform time-harmonic current.


#### Electric Ring Current

The electric ring current with amplitude ``I``, and orientation ``\hat{\bm e}_z`` at position ``\bm r_0`` is assumed to have the current density
```math
\bm{j}_\mathrm{rc} = I \hat{\bm e}_{\varphi'} \delta (r - R_0) \delta (\vartheta - \vartheta_0)\,,
```
where ``\delta`` denotes the Dirac delta distribution and ``R_0 = \sqrt{z_0^2 + a^2}``.

#### Magnetic Ring Current

The magnetic ring current with amplitude ``M``, and orientation ``\hat{\bm e}_z`` at position ``\bm r_0`` is assumed to have the current density
```math
\bm{m}_\mathrm{rc} = M \hat{\bm e}_{\varphi'} \delta (r - R_0) \delta (\vartheta - \vartheta_0)\,.
```

---
## [API](@id rcAPI)

The API provides the following constructors with default values:
```julia
ex = electricRingCurrent(
        embedding  = Medium(ε0, μ0),
        wavenumber = error("missing argument `wavenumber`"),
        amplitude  = 1.0,
        radius     = error("missing argument `radius`"),
        center     = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
        rotation   = SVector{2,typeof(wavenumber)}(0.0, 0.0),
)

ex = magneticRingCurrent(
        embedding  = Medium(ε0, μ0),
        wavenumber = error("missing argument `wavenumber`"),
        amplitude  = 1.0,
        radius     = error("missing argument `radius`"),
        center     = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
        rotation   = SVector{2,typeof(wavenumber)}(0.0, 0.0),
)
```


---
## Radiated Field

The electric field of the electric ring current itself (without scatterer) is computed as a series expansion as defined in [[1, pp. 362ff]](@ref refs).

!!! note
    The fields of the magnetic ring currents are computed via [duality relations](@ref dualityRelations).

!!! note
    It seems there exist other ways of computing the radiated fields. The different formulations should be compared someday.


#### API

The general API is employed:
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
```

---
## Scattered Field

The scattered field computation follows [[1, pp. 368ff]](@ref refs). For the magnetic ring current [duality relations](@ref dualityRelations) are employed. 

!!! note
    Orientation and location of the ring currents are restricted for the scattered field computation!

The ring current has to be oriented such that it is normal to the sphere surface and located outside of the sphere.

!!! warning
    So far the ring current is assumed to be along the ``z``-axis! This is planned to be generalized.

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