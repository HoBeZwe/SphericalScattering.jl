
# Dipoles

```@raw html
<figure>
  <img
    src="assets/Fig_SphereDP.svg"
    alt="Setup"
    width="300" />

  <figcaption>
    Setup for the excitation with a Hertzian or a Fitzgerald dipole.
  </figcaption>
</figure>
<br/>
```

---
## Definition

The dipoles are defined as infinitesimal current elements. Note that the definitions differ from the ones employed in[^1] basically by a factor of ``k``. Hence, they lead to different static fields, i.e., for ``k\rightarrow 0``.

[^1]: Book by Jackson.

#### Hertzian Dipole

The Hertzian dipole with dipole length ``l``, electric current ``I``, and orientation ``\hat{\bm p}`` at position ``\bm r_0`` is assumed to have the current density
```math
\bm{j}_\mathrm{HD} = Il \hat{\bm p} \delta (\bm r - \bm r_0) \,,
```
where ``\delta`` denotes the Dirac delta distribution.

#### Fitzgeral Dipole

The Fitzgerald (magnetic) dipole with dipole length ``l``, magnetic current ``M``, and orientation ``\hat{\bm p}`` at position ``\bm r_0`` is assumed to have the current density
```math
\bm{m}_\mathrm{FD} = Ml \hat{\bm p} \delta (\bm r - \bm r_0) \,.
```

---
## [API](@id dipolesAPI)

The API provides the following constructors with default values:
```julia
ex = HertzianDipole(
        embedding   = Medium(ε0, μ0),
        wavenumber  = error("missing argument `wavenumber`"),
        amplitude   = 1.0,
        center      = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
        orientation = SVector{3,typeof(wavenumber)}(0.0, 0.0, 1.0),
)

ex = FitzgeraldDipole(
        embedding   = Medium(ε0, μ0),
        wavenumber  = error("missing argument `wavenumber`"),
        amplitude   = 1.0,
        center      = SVector{3,typeof(wavenumber)}(0.0, 0.0, 0.0),
        orientation = SVector{3,typeof(wavenumber)}(0.0, 0.0, 1.0),
)
```

!!! warning
    The orientation is so far not automatically normalized!


---
## Radiated Field

The electric field of the Hertzian dipole itself (without scatterer) is
```math
\bm e(\bm r) = Z_\mathrm{F} \cfrac{Il}{4 \pi} \mathrm{e}^{-\mathrm{j} k r}  \left( \cfrac{k}{r}  ((\hat{\bm n} \times \hat{\bm p}) \times \hat{\bm n}) + \left(\cfrac{1}{k r^3} + \cfrac{\mathrm{j}}{r^2} \right)  (3 \hat{\bm n} (\hat{\bm n} \cdot \hat{\bm p}) - \hat{\bm p}) \right)
```
with ``Z_\mathrm{F} = \sqrt{\mu / \varepsilon}`` and the magnetic field
```math
\bm h(\bm r) = \cfrac{Il}{4 \pi} \cfrac{\mathrm{e}^{-\mathrm{j} k r} }{r}  \left(k + \cfrac{1}{\mathrm{j} r} \right) (\hat{\bm n} \times \hat{\bm p})
```
where
```math
\hat{\bm n} = \cfrac{\bm r - \bm r_0}{|\bm r - \bm r_0|}\,.
```
!!! note
    The fields of the Fitzgerald dipole are computed via [duality relations](@ref dualityRelations).


#### API

The general API is employed:
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
```

---
## Scattered Field

!!! note
    Orientation and location of the dipole are restricted for the scattered field computation!

The dipole has to be oriented such that it is normal to the sphere surface and located outside of the sphere.

!!! warning
    So far the dipole is assumed to be along the ``z``-axis! This is planned to be generalized.

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