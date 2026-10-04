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

The scattered field computation follows [bowmanElectromagneticAcousticScattering1970](@cite). It is obtained by a modal series, as for the acoustic [plane wave](@ref ACpwAPI): the addition theorem of the free-space Green's function,
```math
\cfrac{\mathrm{e}^{-\mathrm{j} k |\bm r - \bm r_0|}}{4 \pi |\bm r - \bm r_0|}
    = \cfrac{-\mathrm{j} k}{4 \pi} \sum_{n=0}^\infty (2n+1) j_n(k r_<) h_n^{(2)}(k r_>) P_n(\cos \vartheta) \,,
```
with ``r_< = \min(r, r_0)``, ``r_> = \max(r, r_0)``, ``r_0 = |\bm r_0|``, the Legendre polynomials ``P_n`` and ``\vartheta`` measured from the direction ``\hat{\bm r}_0`` towards the monopole, expands the incident pressure in the spherical waves about the center of the sphere. Since the monopole lies outside the sphere, ``r_> = r_0`` holds on its surface, where the boundary condition is imposed. The scattered pressure therefore reads
```math
p^\mathrm{sc}(\bm r) = \cfrac{-\mathrm{j} k a}{4 \pi} \sum_{n=0}^\infty (2n+1) h_n^{(2)}(k r_0) \, b_n \, h_n^{(2)}(k r) P_n(\cos \vartheta)
```
with the very same scattering coefficients
```math
b_n = -\cfrac{j_n^\prime(k r_\mathrm{s})}{h_n^{(2)\prime}(k r_\mathrm{s})}
\qquad \text{and} \qquad
b_n = -\cfrac{j_n(k r_\mathrm{s})}{h_n^{(2)}(k r_\mathrm{s})}
```
for a sound-hard and a sound-soft sphere of radius ``r_\mathrm{s}`` as for the plane wave [bowmanElectromagneticAcousticScattering1970](@cite). Applying the boundary condition term by term leaves the radial structure of the series untouched, hence the scattering coefficients do not depend on the excitation; only the coefficients of the incident expansion do.

!!! note
    In contrast to the electromagnetic case the series starts at ``n=0``: the monopole term contributes and dominates the low-frequency limit.

!!! note
    No rotation of the coordinate system is required for an arbitrary position of the monopole. The scattered pressure is a scalar and rotationally symmetric about ``\hat{\bm r}_0``, so that it depends on the observation point only via ``r`` and ``\cos \vartheta = \hat{\bm r}_0 \cdot \hat{\bm r}``.

!!! warning
    The monopole has to be located outside the sphere, as the expansion above assumes. This is checked during the computation, an error being thrown otherwise.

!!! tip
    The series converges like ``(r_\mathrm{s}^2 / (r_0 r))^n``. A monopole very close to the surface therefore requires many terms, which may exceed the numerical range of the spherical Hankel functions before the series converges; a message is printed should this happen.

#### API

The general API is employed:
```julia
p  = scatteredfield(sp, ex, Pressure(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))

γ₀ = scatteredfield(sp, ex, PressureTrace(point_cart))

γ₁ = scatteredfield(sp, ex, PressureNormalGradient(point_cart))
```
where `sp` is a [`HardSphere`](@ref) or a [`SoftSphere`](@ref). A [`Spheroid`](@ref) or a [`Disc`](@ref) is scattered from as well, by a series in the oblate spheroidal wave functions instead; see [Spheroid and Disc](@ref ACspheroidAPI).

!!! note
    The traces are evaluated on the surface of the sphere: only the direction of each location is taken into account, the radial coordinate being replaced by the radius of the sphere. Hence, the locations may also be given by the points of a faceted surface mesh, which do not lie exactly on the sphere. Normals other than ``\hat{\bm r}`` are supported, the tangential part of the gradient being included.


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