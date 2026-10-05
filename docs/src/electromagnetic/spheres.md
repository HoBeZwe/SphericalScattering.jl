
# Spheres

!!! note
    In all of the following setups the sphere is embedded in a homogeneous background medium with permeability ``\mu`` and permittivity ``\varepsilon``.
    This medium is defined only by the chosen excitation as a [`Medium(ε, μ)`](@ref). 


---
## PEC/PMC Sphere

The perfectly electrically conducting (PEC) or perfectly magnetically conducting (PMC) sphere has radius ``r`` and is assumed to be located in the origin. It is defined by [`PECSphere`](@ref).
```@raw html
<div align="center">
<img src="../../assets/PECsphere.svg" width="300"/>
</div>
<br/>
```

#### [API](@id pecAPI)

```@docs
PECSphere
```
`PECSphere` is an alias of `Sphere{PEC}`: a sphere carrying the boundary condition [`PEC`](@ref) as its type parameter, see [`Scatterer`](@ref).


---
## Dielectric Sphere

The dielectric sphere has radius ``r`` and is assumed to be located in the origin. It is defined by [`DielectricSphere`](@ref), where the filling [`Medium(εᵢ, μᵢ)`](@ref) with permeability ``\mu_\mathrm{i}`` and permittivity ``\varepsilon_\mathrm{i}`` has to be defined. 
```@raw html
<div align="center">
<img src="../../assets/DielectricSphere.svg" width="300"/>
</div>
<br/>
```

#### [API](@id dielecAPI)

```@docs
DielectricSphere
```
Here `radius` is a Float and `filling` is of type [`Medium(εᵢ, μᵢ)`](@ref).

`DielectricSphere` is an alias of `Sphere{<:Dielectric}`: the filling is carried by the boundary condition [`Dielectric`](@ref), which, in contrast to [`PEC`](@ref), requires data and is therefore stored as a value as well, `sp.boundary.filling`.


---
## Layered Dielectric Sphere

The layered dielectric sphere has radii ``[r_1, r_2, \dots, r_N]`` and is assumed to be located in the origin. It is defined by [`LayeredSphere`](@ref), where the vector of fillings [[`Medium(ε₁, μ₁)`](@ref), [`Medium(ε₂, μ₂)`](@ref), ..., [`Medium(εN, μN)`](@ref)] with permeability ``\mu_n`` and permittivity ``\varepsilon_n`` has to be defined.
```@raw html
<div align="center">
<img src="../../assets/LayeredSphere.svg" width="300"/>
</div>
<br/>
```

#### [API](@id mlDielecAPI)

```@docs
LayeredSphere
```
with, e.g., `radii = SVector(0.25, 0.5, 1.0)` and `filling = SVector(Medium(ε1, μ1), Medium(ε2, μ2), Medium(ε3, μ3))`.

`LayeredSphere` is an alias of `Sphere{<:Layered{<:Dielectric}}`. Seen from outside, the shells form a condition on the outer surface: the outermost radius ``r_N`` is the radius of the [`Sphere`](@ref SphericalScattering.Sphere), while the inner interfaces, the fillings of the shells, and the dielectric core are carried by the boundary condition [`Layered`](@ref).


---
## Layered Dielectric Sphere with PEC Core

The layered dielectric sphere has radii ``[r_1, r_2, \dots, r_{N+1}]`` and is assumed to be located in the origin. It is defined by [`LayeredSpherePEC`](@ref), where the vector of fillings [[`Medium(ε₁, μ₁)`](@ref), [`Medium(ε₂, μ₂)`](@ref), ..., [`Medium(εN, μN)`](@ref)] with permeability ``\mu_n`` and permittivity ``\varepsilon_n`` has to be defined.
```@raw html
<div align="center">
<img src="../../assets/LayeredSpherePEC.svg" width="300"/>
</div>
```

#### [API](@id mlDielecPecAPI)

```@docs
LayeredSpherePEC
```
with, e.g., `radii = SVector(0.25, 0.5, 1.0)` and `filling = SVector(Medium(ε1, μ1), Medium(ε2, μ2))`.

`LayeredSpherePEC` is an alias of `Sphere{<:Layered{PEC}}`: as for the [`LayeredSphere`](@ref), the outermost radius is that of the sphere, and the condition [`Layered`](@ref) carries the rest, here with a [`PEC`](@ref) core.


---
## [Dielectric Sphere with Thin Impedance Layer](@id dielecimped)

The dielectric sphere with a thin impedance layer of thickness ``t`` has radius ``r`` and is assumed to be located in the origin. It is defined by [`DielectricSphereThinImpedanceLayer`](@ref). Unlike the LayeredSphere model, the solution is obtained by using an approximation: it is assumed that the impedance is so high that the displacement field is purely radial (see [jonesElectromechanicsParticles1995; pp. 230ff](@cite)). This leads to a potential drop across the thin layer, while the displacement field is constant in radial direction. In addition to the filling [`Medium(εᵢ, μᵢ)`](@ref), the impedance layer must be specified, both the [`Medium(εₜ, μₜ)`](@ref) and its `thickness`.  
```@raw html
<div align="center">
<img src="../../assets/ImpedanceLayer.svg" width="300"/>
</div>
<br/>
```

!!! note
    This configuration is (at least so far) only intended for the [uniform static field](@ref uniformEx) excitation.

#### [API](@id dielecimpedAPI)

```@docs
DielectricSphereThinImpedanceLayer
```
Here `radius` and `thickness` are a Floats, `filling` and `thinlayer` are of type [`Medium`](@ref).

`DielectricSphereThinImpedanceLayer` is an alias of `Sphere{<:ThinImpedanceLayer}`: the coating, being thin, is modelled as an effective condition on the surface of the sphere, see [`ThinImpedanceLayer`](@ref), which carries the `thickness`, the `thinlayer`, and the `filling`.
