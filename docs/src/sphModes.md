
# Spherical Modes

```@raw html
<figure>
  <img
    src="../assets/PECsphere.svg"
    alt="Setup"
    width="300" />

  <figcaption>
    Setup for the excitation with spherical modes. TODO: more suitable image.
  </figcaption>
</figure>
<br/>
```

!!! info
    So far only single modes (not their superposition) are implemented.


---
## Definition

The spherical modes are defined following the conventions of [[5]](@ref refs), however adapted to the time convention ``\mathrm{e}^{\,\mathrm{j}\omega t}``.

#### TE Modes

The functions
```math
\begin{aligned}
\bm{f}_{1mn}^{(c)}(r,\vartheta, \varphi) = \cfrac{1}{2\pi} \cfrac{1}{n(n+1)}{\left(-\frac{m}{|m|}\right)}^m 
													&\left(\mathrm{h}_n^{(c)}(kr)\frac{\mathrm{j}m\bar{\mathrm{P}}_n^{|m|}(\cos\vartheta)}{\sin\vartheta}\mathrm{e}^{\mathrm{j}m\varphi} \bm{e}_\vartheta   	\right.	 \\	
													&~~\left. -\mathrm{h}_n^{(c)}(kr)\cfrac{\mathrm{d} \bar{\mathrm{P}}_n^{|m|}(\cos\vartheta) }{\mathrm{d}\vartheta} \mathrm{e}^{\mathrm{j}m\varphi} \bm{e}_\varphi  	~\right)
\end{aligned}
```
form the ``\mathrm{TE}_{mn}`` modes, where ``\mathrm{P}_n^m`` denote the associated Legendre polynomials and ``\mathrm{h}_n^{(c)}(x)`` denote the spherical Hankel functions.


#### TM Modes

The functions
```math
\begin{aligned}
\bm{f}_{2mn}^{(c)}(r,\vartheta, \varphi) = \cfrac{1}{2\pi} \cfrac{1}{n(n+1)}{\left(-\frac{m}{|m|}\right)}^m \Bigg(&\cfrac{n(n+1)}{kr}\mathrm{h}_n^{(c)}   \bar{\mathrm{P}}_n^{|m|}(\cos\vartheta)	\mathrm{e}^{\mathrm{j}m\varphi} \bm{e}_r \Bigg. \\	
					&+ \cfrac{1}{kr} \cfrac{\mathrm{d}\left( kr \mathrm{h}_n^{(c)}(kr)\right)}{\mathrm{d}(kr)} \cfrac{\mathrm{d} \bar{\mathrm{P}}_n^{|m|}(\cos\vartheta) }{\mathrm{d}\vartheta} \mathrm{e}^{\mathrm{j}m\varphi} \bm{e}_\vartheta  \\
					&\Bigg. \frac{1}{kr} \cfrac{\mathrm{d}\left( kr \mathrm{h}_n^{(c)}(kr)\right)}{\mathrm{d}(kr)}\cfrac{\mathrm{j}m\bar{\mathrm{P}}_n^{|m|}(\cos\vartheta)}{\sin\vartheta}\mathrm{e}^{\mathrm{j}m\varphi} \bm{e}_\varphi	\Bigg) \,\text{.}
\end{aligned}
```
form the ``\mathrm{TM}_{mn}`` modes.

!!! note
    By ``c=1`` inward travelling waves and by ``c=2`` outward travelling waves are denoted.
    
    The change of time convention from [[5]](@ref refs) to the one of this package is achieved solely by interchanging these values. 

---
## [API](@id modesAPI)

The API provides the following constructors with default values:
```@docs
SphericalModeTE
SphericalModeTM
```

---
## Incident Field

The electric field is given by
```math
\bm e^\mathrm{TE,inc}_{mn} = k \sqrt{Z_\mathrm{F}} \bm{f}_{1mn}^{(1)}
```
and 
```math
\bm e^\mathrm{TM,inc}_{mn} = k \sqrt{Z_\mathrm{F}} \bm{f}_{2mn}^{(1)}
```

The magnetic field is given by
```math
\bm e^\mathrm{TE,inc}_{mn} = \cfrac{\mathrm{j} k}{\sqrt{Z_\mathrm{F}}} \bm{f}_{2mn}^{(1)}
```
and 
```math
\bm e^\mathrm{TM,inc}_{mn} = \cfrac{\mathrm{j} k}{\sqrt{Z_\mathrm{F}}} \bm{f}_{1mn}^{(1)}
```


!!! note
    The modes are assumed to travel towards the origin.



#### API

The general API is employed:
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
```

!!! note
    The far-field for a spherical mode makes only sense when travelling outwards.


---
## Scattered Field

Matching incoming and outcoming waves to fulfull the boundary condition ``\bm e_\mathrm{tan} = \bm 0`` yields the scattering coefficients
```math
\xi_\mathrm{TE} = -\cfrac{\mathrm{H}^{(2)}_{n + 0.5}(k a)}{\mathrm{H}^{(1)}_{n + 0.5}(k a)}
```
and
```math
\xi_\mathrm{TM} = -\cfrac{\mathrm{H}'^{(2)}_{n + 0.5}(k a)}{\mathrm{H}'^{(1)}_{n + 0.5}(k a)}
```
where ``\mathrm{H}^{(\nu)}_{n}(x)`` denotes the Hankel function of ``\nu``-th kind and ``n``-th order.

The scattered fields ``\bm e^\mathrm{sc}`` are then given by
```math
\bm e_\mathrm{TE}^\mathrm{sc} = \xi_\mathrm{TE} k \sqrt{Z_\mathrm{F}} \bm{f}_{1mn}^{(2)}
```
and
```math
\bm e_\mathrm{TM}^\mathrm{sc} = \xi_\mathrm{TM} k \sqrt{Z_\mathrm{F}} \bm{f}_{2mn}^{(2)} \,.
```

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
    The total far-field for a spherical mode excitation is not defined.
