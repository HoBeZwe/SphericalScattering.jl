# Monopole

---
## Definition

The acoustic monopole (point source) with amplitude ``a`` at position ``\bm r_0`` is defined by the pressure field [bowmanElectromagneticAcousticScattering1970](@cite)
```math
p_\mathrm{M}(\bm r) = \cfrac{a}{4 \pi} \cfrac{\mathrm{e}^{-\mathrm{j} k R}}{R} \,, \qquad R = |\bm r - \bm r_0| \,,
```
that is, by the free-space Green's function of the Helmholtz equation scaled by the amplitude. The field is singular at the position of the monopole.

!!! note
    The background medium is defined by the chosen excitation as a [`Medium(ε, μ)`](@ref), just as for the electromagnetic excitations. For acoustics the two parameters are to be read via the analogy ``\mu \rightarrow \rho`` (mass density) and ``\varepsilon \rightarrow 1 / (\rho c^2)`` (compressibility), so that the speed of sound is ``c = 1 / \sqrt{\varepsilon \mu}`` and the wave impedance is ``\rho c = \sqrt{\mu / \varepsilon}``. Since the default is free space, a medium matching the fluid of interest should be provided.

!!! note
    Throughout this package the time convention ``\mathrm{e}^{\mathrm{j} \omega t}`` is employed, hence the spherical Hankel functions of the second kind ``h_n^{(2)}`` describe outgoing waves.


---
## [API](@id ACpointAPI)

The API provides the following constructor with default values:
```@docs
Acoustic.monopole
```

The constructors of the acoustic excitations are collected in the [`Acoustic` submodule](@ref ACsubmodule), so that their names do not clash with the electromagnetic ones.


---
## Radiated Field

The pressure of the monopole itself (without scatterer) is as given above. Since it depends on the position solely via the distance ``R`` to the monopole, its gradient is radial with respect to the monopole,
```math
\nabla p_\mathrm{M} = \cfrac{\partial p_\mathrm{M}}{\partial R} \hat{\bm R}
                    = -\cfrac{a}{4 \pi} (1 + \mathrm{j} k R) \cfrac{\mathrm{e}^{-\mathrm{j} k R}}{R^2} \hat{\bm R} \,,
\qquad \hat{\bm R} = \cfrac{\bm r - \bm r_0}{|\bm r - \bm r_0|} \,,
```
so that the two Cauchy data (the surface traces) on a surface with outward unit normal ``\hat{\bm n}`` read
```math
\gamma_0 p_\mathrm{M} = p_\mathrm{M}
\qquad \text{and} \qquad
\gamma_1 p_\mathrm{M} = \hat{\bm n} \cdot \nabla p_\mathrm{M}
                      = -\cfrac{a}{4 \pi} (1 + \mathrm{j} k R) \cfrac{\mathrm{e}^{-\mathrm{j} k R}}{R^2} (\hat{\bm R} \cdot \hat{\bm n}) \,.
```
For ``k R \gg 1`` the Neumann trace approaches ``-\mathrm{j} k (\hat{\bm R} \cdot \hat{\bm n}) p_\mathrm{M}``, the local plane-wave result.

Unlike a plane wave, a monopole does possess a far field. With ``|\bm r - \bm r_0| = r - \hat{\bm r} \cdot \bm r_0 + \mathcal{O}(1/r)`` it is
```math
p_{\mathrm{M},\infty}(\hat{\bm r}) = \lim_{r \rightarrow \infty} r \, \mathrm{e}^{\mathrm{j} k r} p_\mathrm{M}(\bm r)
                                   = \cfrac{a}{4 \pi} \mathrm{e}^{\mathrm{j} k \hat{\bm r} \cdot \bm r_0} \,,
```
that is, the factor ``\mathrm{e}^{-\mathrm{j} k r} / r`` is omitted, as for the electromagnetic excitations. Its magnitude is the same in all directions, as it has to be for a point source: the position of the monopole enters the phase alone.

#### API

The general API is employed:
```julia
p  = field(ex, Pressure(point_cart))

FF = field(ex, FarField(point_cart))

γ₀ = field(ex, PressureTrace(point_cart))

γ₁ = field(ex, PressureNormalGradient(point_cart))
```
Without further specification the outward normal ``\hat{\bm n} = \hat{\bm r}`` of a sphere centered in the origin is assumed at every location. One normal vector per location may be provided instead, which is required, e.g., for the facets of a surface mesh:
```julia
γ₁ = field(ex, PressureNormalGradient(point_cart, normals))
```
See [`PressureNormalGradient`](@ref) for details.


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

!!! warning
    The monopole has to be located outside the scatterer, as the expansion of its field assumes. This is checked during the computation, an error being thrown otherwise.


---
## Total Field

#### API

The general API is employed:
```julia
p  = field(sp, ex, Pressure(point_cart))

FF = field(sp, ex, FarField(point_cart))

γ₀ = field(sp, ex, PressureTrace(point_cart))

γ₁ = field(sp, ex, PressureNormalGradient(point_cart))
```

!!! note
    In contrast to the plane wave, the total far-field is defined, since the incident far-field is.

!!! tip
    The total field is the Green's function of the exterior problem and hence symmetric in the position of the monopole and the observation point. Together with the boundary conditions — ``\gamma_1 p`` vanishing on a sound-hard and ``\gamma_0 p`` on a sound-soft sphere — this is a convenient check of the implementation.