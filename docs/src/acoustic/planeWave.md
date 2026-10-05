# Plane Wave

---
## Definition

An acoustic plane wave with amplitude ``a`` and direction of propagation ``\hat{\bm d}`` (vectors with a hat denote unit vectors) is defined by the pressure field [bowmanElectromagneticAcousticScattering1970](@cite)
```math
p_\mathrm{PW}(\bm r) = a \, \mathrm{e}^{-\mathrm{j} k \hat{\bm d} \cdot \bm r}
```
with the wavenumber ``k = \omega \sqrt{\varepsilon \mu}``.

!!! note
    The background medium is defined by the chosen excitation as a [`Medium(ε, μ)`](@ref), just as for the electromagnetic excitations. For acoustics the two parameters are to be read via the analogy ``\mu \rightarrow \rho`` (mass density) and ``\varepsilon \rightarrow 1 / (\rho c^2)`` (compressibility), so that the speed of sound is ``c = 1 / \sqrt{\varepsilon \mu}`` and the wave impedance is ``\rho c = \sqrt{\mu / \varepsilon}``. Since the default is free space, a medium matching the fluid of interest should be provided.

!!! note
    Throughout this package the time convention ``\mathrm{e}^{\mathrm{j} \omega t}`` is employed, hence the spherical Hankel functions of the second kind ``h_n^{(2)}`` describe outgoing waves.


---
## [API](@id ACpwAPI)

The API provides the following constructor with default values:
```@docs
Acoustic.planeWave
```

!!! tip
    The `direction` vector is automatically normalized to a unit vector during the initialization.

#### [The `Acoustic` Submodule](@id ACsubmodule)

The constructors of the acoustic excitations are collected in a submodule, so that their names do not clash with the electromagnetic ones:
```@docs
Acoustic
```


---
## Incident Field

The pressure of the plane wave is as given above. On a surface with outward unit normal ``\hat{\bm n}`` the two Cauchy data (the surface traces) are the Dirichlet trace
```math
\gamma_0 p_\mathrm{PW} = p_\mathrm{PW}
```
and, since ``\nabla p_\mathrm{PW} = -\mathrm{j} k \hat{\bm d} \, p_\mathrm{PW}``, the Neumann trace
```math
\gamma_1 p_\mathrm{PW} = \hat{\bm n} \cdot \nabla p_\mathrm{PW} = -\mathrm{j} k (\hat{\bm d} \cdot \hat{\bm n}) \, p_\mathrm{PW} \,.
```

#### API

The general API is employed:
```julia
p  = field(ex, Pressure(point_cart))

γ₀ = field(ex, PressureTrace(point_cart))

γ₁ = field(ex, PressureNormalGradient(point_cart))
```
Without further specification the outward normal ``\hat{\bm n} = \hat{\bm r}`` of a sphere centered in the origin is assumed at every location. One normal vector per location may be provided instead, which is required, e.g., for the facets of a surface mesh:
```julia
γ₁ = field(ex, PressureNormalGradient(point_cart, normals))
```
See [`PressureNormalGradient`](@ref) for details.

!!! note
    The far-field of a plane wave is not defined.


---
## Scattered Field

The scattered field computation follows [bowmanElectromagneticAcousticScattering1970](@cite). For the sound-hard and the sound-soft sphere it is obtained by a series in the spherical wave functions, see [Acoustic Spheres](@ref acSphereSeries); for the spheroids and the disc by a series in the spheroidal wave functions, see [Spheroidal Solution](@ref ACspheroidSolution).

#### API

The general API is employed:
```julia
p  = scatteredfield(sp, ex, Pressure(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))

γ₀ = scatteredfield(sp, ex, PressureTrace(point_cart))

γ₁ = scatteredfield(sp, ex, PressureNormalGradient(point_cart))
```
where `sp` is a [`HardSphere`](@ref), a [`SoftSphere`](@ref), a [`Spheroid`](@ref), or a [`Disc`](@ref). The traces and the far field are described under [Quantities](@ref quantitiesConcept).

!!! warning
    For a spheroid, the outward normals have to be provided, `PressureNormalGradient(sp, point_cart)`, see [Surface Traces](@ref).

!!! note
    The pressure vanishes inside the scatterer, whereas the far field is determined by the direction of observation alone and is, therefore, not suppressed for locations inside the scatterer.


---
## Total Field

#### API

The general API is employed:
```julia
p  = field(sp, ex, Pressure(point_cart))

γ₀ = field(sp, ex, PressureTrace(point_cart))

γ₁ = field(sp, ex, PressureNormalGradient(point_cart))
```

!!! note
    The total far-field is not defined (since the incident far-field is not defined).

!!! tip
    The total traces are a convenient check of the implementation: ``\gamma_1 p`` vanishes on a sound-hard sphere and ``\gamma_0 p`` vanishes on a sound-soft sphere.
