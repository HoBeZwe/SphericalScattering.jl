
# Manual

The basic building blocks are introduced in the following simple example; more details are provided afterwards:

---
## Introductory Example: Plane Wave



```julia
using SphericalScattering, StaticArrays

# define excitation: plane wave travelling in negative z-direction with x-polarization
ex = planeWave(wavenumber=30.0) # ≈ 10 MHz

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

### Ring currents
```julia
ex = electricRingCurrent(
        embedding   = Medium(ε, μ),
        wavenumber  = 30.0,
        amplitude   = 1.0,
        radius      = 0.5,
        center      = SVector(0.0,0.0,0.0),
        rotation    = SVector(0.0,0.0)
)

ex = magneticRingCurrent(
        embedding   = Medium(ε, μ),
        wavenumber  = 30.0,
        amplitude   = 1.0,
        radius      = 0.5,
        center      = SVector(0.0,0.0,0.0),
        rotation    = SVector(0.0,0.0)
)
```

### Spherical modes

```julia
ex = SphericalModeTE(
        embedding   = Medium(ε, μ),
        wavenumber  = 30.0,
        amplitude   = 1.0,
        m           = 0,
        n           = 1,
        c           = 1,
        center      = SVector(0.0,0.0,0.0),
        orientation = SVector(0.0,0.0,1.0)
)

ex = SphericalModeTM(
        embedding   = Medium(ε, μ),
        wavenumber  = 30.0,
        amplitude   = 1.0,
        m           = 0,
        n           = 1,
        c           = 1,
        center      = SVector(0.0,0.0,0.0),
        orientation = SVector(0.0,0.0,1.0)
)
```

### Uniform Static Field

```julia
ex = UniformField(
        embedding    = Medium(ε0, μ0),
        amplitude    = 1.0,
        direction    = SVector(1.0,0.0,0.0)
)
```

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

#### Incident Fields

For each excitation the far-field, the electric, and the magnetic near-field without a scatterer can be determined: 
```julia
E  = field(ex, ElectricField(point_cart))

H  = field(ex, MagneticField(point_cart))

FF = field(ex, FarField(point_cart))
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

E = scatteredfield(sp, ex, ElectricField(point_cart))
```

---
## Conversion Between Basis

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