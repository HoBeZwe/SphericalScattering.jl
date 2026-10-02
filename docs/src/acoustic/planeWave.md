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

The scattered field computation follows [bowmanElectromagneticAcousticScattering1970](@cite). It is obtained by a modal series, as for the electromagnetic excitations: expanding the incident pressure in the spherical waves about the center of the sphere,
```math
p_\mathrm{PW}(\bm r) = a \sum_{n=0}^\infty (2n+1) (-\mathrm{j})^n j_n(k r) P_n(\cos \vartheta) \,,
```
where ``\vartheta`` is measured from the direction of incidence ``\hat{\bm d}`` and ``P_n`` denotes the Legendre polynomials, the scattered pressure follows as
```math
p^\mathrm{sc}(\bm r) = a \sum_{n=0}^\infty (2n+1) (-\mathrm{j})^n b_n h_n^{(2)}(k r) P_n(\cos \vartheta) \,.
```
Applying the boundary condition term by term yields the scattering coefficients
```math
b_n = -\cfrac{j_n^\prime(k r_\mathrm{s})}{h_n^{(2)\prime}(k r_\mathrm{s})}
\qquad \text{and} \qquad
b_n = -\cfrac{j_n(k r_\mathrm{s})}{h_n^{(2)}(k r_\mathrm{s})}
```
for a sound-hard and a sound-soft sphere of radius ``r_\mathrm{s}``, respectively: on a sound-hard surface the normal velocity, and hence the radial derivative of the total pressure, vanishes, whereas on a sound-soft (pressure release) surface the total pressure vanishes [bowmanElectromagneticAcousticScattering1970](@cite).

!!! note
    In contrast to the electromagnetic case the series starts at ``n=0``: the monopole term contributes and dominates the low-frequency limit.

!!! note
    No rotation of the coordinate system is required for an arbitrary direction of incidence. The scattered pressure is a scalar and rotationally symmetric about ``\hat{\bm d}``, so that it depends on the observation point only via ``r`` and ``\cos \vartheta = \hat{\bm d} \cdot \hat{\bm r}``.

The far field is defined as
```math
p^\mathrm{sc}_\infty(\hat{\bm r}) = \lim_{r \rightarrow \infty} r \, \mathrm{e}^{\mathrm{j} k r} p^\mathrm{sc}(\bm r)
                                  = \cfrac{\mathrm{j} a}{k} \sum_{n=0}^\infty (2n+1) b_n P_n(\cos \vartheta) \,,
```
that is, the factor ``\mathrm{e}^{-\mathrm{j} k r} / r`` is omitted, as for the electromagnetic excitations. Note that the order-dependent factors cancel, since ``(-\mathrm{j})^n \mathrm{j}^{n+1} = \mathrm{j}`` holds for every ``n``.

#### API

The general API is employed:
```julia
p  = scatteredfield(sp, ex, Pressure(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))

γ₀ = scatteredfield(sp, ex, PressureTrace(point_cart))

γ₁ = scatteredfield(sp, ex, PressureNormalGradient(point_cart))
```
where `sp` is a [`HardSphere`](@ref) or a [`SoftSphere`](@ref).

!!! note
    The traces are evaluated on the surface of the sphere: only the direction of each location is taken into account, the radial coordinate being replaced by the radius of the sphere. Hence, the locations may also be given by the points of a faceted surface mesh, which do not lie exactly on the sphere.

!!! tip
    For a normal ``\hat{\bm n} \neq \hat{\bm r}`` the Neumann trace of the scattered field picks up the tangential part of the gradient as well:
    ```math
    \hat{\bm n} \cdot \nabla p^\mathrm{sc} = (\hat{\bm n} \cdot \hat{\bm r}) \cfrac{\partial p^\mathrm{sc}}{\partial r}
        + (\hat{\bm n} \cdot \hat{\bm \vartheta}) \cfrac{1}{r} \cfrac{\partial p^\mathrm{sc}}{\partial \vartheta} \,.
    ```
    Both contributions are included, so that providing normals is supported here as well.

!!! note
    The pressure vanishes inside the sphere, whereas the far field is determined by the direction of observation alone and is, therefore, not suppressed for locations inside the sphere.


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
