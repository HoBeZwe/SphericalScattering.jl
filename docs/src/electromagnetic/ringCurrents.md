
# Ring Currents

```@raw html
<figure>
  <img
    src="../assets/Fig_SphereRC.svg"
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

The electric ring current with amplitude ``I``, and orientation ``\hat{\bm p} = \hat{\bm e}_z`` with its center at position ``\bm r_0 = (0,0,z_0)`` is assumed to have the current density
```math
\bm{j}_\mathrm{rc} = I \hat{\bm e}_{\varphi'} \delta (r - R_0) \delta (\vartheta - \vartheta_0)\,,
```
where ``\delta`` denotes the Dirac delta distribution and ``R_0 = \sqrt{z_0^2 + a^2}``.

#### Magnetic Ring Current

The magnetic ring current with amplitude ``M``, and orientation ``\hat{\bm p} = \hat{\bm e}_z`` at position ``\bm r_0 = (0,0,z_0)`` is assumed to have the current density
```math
\bm{m}_\mathrm{rc} = M \hat{\bm e}_{\varphi'} \delta (r - R_0) \delta (\vartheta - \vartheta_0)\,.
```

!!! note
    Arbitrary center positions ``\bm r_0`` and orientations ``\hat{\bm p}`` are computed by [rotations](@ref rotationDetails) and translations of the corresponding fields.



---
## [API](@id rcAPI)

The API provides the following constructors with default values:
```@docs
electricRingCurrent
magneticRingCurrent
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

!!! tip
    For the scattered field computation the orientation of the ring current is restricted: the orientation vector has to be perpendicular to the surface of the sphere. 

    This can be achieved, e.g., by defining the orientation first
    ```julia
    orientation = normalize(SVector(0.0,1.0,1.0))
    ```
    and then setting the center position as a multiple of the orientaion:
    ```julia
    center = 2.0 * orientation
    ```

!!! note
    Internal details of the computations: Following [[1, pp. 368ff]](@ref refs) the orientation vectors of the ring currents are initially assumed to be aligned with the ``z``-axis. Arbitrary centers and orientations (forming a valid pair) are obtained via [rotations](@ref rotationDetails).

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

FF = field(sp, ex, FarField(point_cart))
```