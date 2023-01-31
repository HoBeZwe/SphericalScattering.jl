
# Manual

The basic building blocks are introduced in the following simple example; more details are provided afterwards:

---
## Introductory Example: Plane Wave



```julia
using SphericalScattering, StaticArrays

# define excitation: plane wave travelling in positive z-direction with x-polarization
ex = planeWave(frequency=10e6) # Hz

# define scatterer: PEC sphere
sp = PECSphere(radius = 1.0)

# define an observation point
point_cart = [SVector(2.0, 2.0, 3.2)] 

# compute scattered fields
E  = scatteredfield(sp, ex, ElectricField(point_cart))
H  = scatteredfield(sp, ex, MagneticField(point_cart))
FF = scatteredfield(sp, ex, FarField(point_cart))
```

---
## Defining Observation Points

In order to define points the [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) package has to be used.
```julia
using StaticArrays

# defining a single point
point_cart = [SVector(2.0, 2.0, 3.2)] 

# defining multiple points (along a line)
point_cart = [SVector(5.0, 5.0, z) for z in -2:0.2:2]
```

!!! info
    A Cartesian basis is used for all coordinates and field components!


---
## Defining an Excitation

For all available excitations a simple constructor with keyword arguments and default values is available. For more details see the APIs of the 

- [Plane wave](@ref pwAPI)
- [Dipoles](@ref dipolesAPI)
- [Ring currents](@ref rcAPI)
- [Spherical modes](@ref modesAPI)
- [Uniform static field](@ref uniformAPI)


---
## Defining a Scatterer

For all available scatteres a simple constructor with keyword arguments and default values is available. For more details see the APIs of the 

- [PEC sphere](@ref pecAPI)
- [PMC sphere](@ref pecAPI)
- [Dielectric sphere](@ref dielecAPI)
- [Multilayer dielectric sphere](@ref mlDielecAPI)
- [Multilayer dielectric sphere with PEC core](@ref mlDielecPecAPI)


---
## Computing Fields

The incident, scattered, and total fields can be computed.
```@raw html
<br/>
```

#### Incident Fields

For each excitation the far-field, the electric, and the magnetic near-field without a scatterer can be determined: 
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
```

For the uniform field excitation, only the electric field as well as the scalar potential can be calculated:
```julia
Φ = field(ex, ScalarPotential(point_cart))
```
```@raw html
<br/>
```

#### Scattered Fields

For each excitation the scattered fields from a given sphere can be determined as:

```julia
E  = scatteredfield(sp, ex, ElectricField(point_cart))

H  = scatteredfield(sp, ex, MagneticField(point_cart))

FF = scatteredfield(sp, ex, FarField(point_cart))

```
For the uniform field excitation, only the electric field as well as the scalar potential can be calculated:
```julia
Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))
```
```@raw html
<br/>
```

#### Total Fields

For each excitation the total fields in the presence of a given sphere can be determined as:

```julia
E  = field(sp, ex, ElectricField(point_cart))

H  = field(sp, ex, MagneticField(point_cart))

FF = field(sp, ex, FarField(point_cart))

```
For the uniform field excitation, only the electric field as well as the scalar potential can be calculated:
```julia
Φ = field(sp, ex, ScalarPotential(point_cart))
```


---
## Conversion Between Bases

Methods are provided to convert between Cartesian and spherical coordinates:

```julia
point_cart = SphericalScattering.sph2cart.(point_sph)

point_sph  = SphericalScattering.cart2sph.(point_cart)
```

Converting fields:

```julia
F_cart = SphericalScattering.convertSpherical2Cartesian.(F_sph,  point_sph)

F_sph  = SphericalScattering.convertCartesian2Spherical.(F_cart, point_sph)
```