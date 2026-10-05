
# General Usage

Every computation combines three building blocks: an excitation, a scatterer, and a quantity to be evaluated at given locations. They are introduced by two simple examples, one for each physics; more details are provided afterwards.


---
## Introductory Examples

#### Electromagnetic: Plane Wave and PEC Sphere

```@example introductory
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
nothing # hide
```

#### Acoustic: Plane Wave and Sound-Hard Sphere

```@example introductory
# define the medium: air, with compressibility 1 / (ρc²) and mass density ρ
air = Medium(1 / (1.2 * 343.0^2), 1.2)

# define excitation: plane wave travelling in positive z-direction
ex = Acoustic.planeWave(frequency=1e3, embedding=air) # Hz

# define scatterer: sound-hard sphere
sp = HardSphere(radius = 0.1)

# define an observation point
point_cart = [SVector(0.2, 0.2, 0.3)]

# compute scattered pressure
p  = scatteredfield(sp, ex, Pressure(point_cart))
FF = scatteredfield(sp, ex, FarField(point_cart))
nothing # hide
```

!!! note
    The background medium is defined by the excitation as a [`Medium(ε, μ)`](@ref), for both physics. For acoustics, the two parameters are read as the compressibility ``\varepsilon \rightarrow 1 / (\rho c^2)`` and the mass density ``\mu \rightarrow \rho``, see the [acoustic plane wave](@ref ACpwAPI). Since the default is free space, a medium matching the fluid of interest should be provided.


---
## Defining Observation Points

In order to define points the [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) package has to be used.
```@example introductory
using StaticArrays

# defining a single point
point_cart = [SVector(2.0, 2.0, 3.2)] 

# defining multiple points (along a line)
point_cart = [SVector(5.0, 5.0, z) for z in -2:0.2:2]
nothing # hide
```

!!! info
    A Cartesian basis is used for all coordinates and field components!


---
## Defining an Excitation

For all available excitations a simple constructor with keyword arguments and default values is available:

| physics | excitation | constructor | details |
|:------- |:---------- |:----------- |:------- |
| electromagnetic | plane wave | `planeWave` | [plane wave](@ref pwAPI) |
| electromagnetic | Hertzian and Fitzgerald dipole | `HertzianDipole`, `FitzgeraldDipole` | [dipoles](@ref dipolesAPI) |
| electromagnetic | electric and magnetic ring current | `electricRingCurrent`, `magneticRingCurrent` | [ring currents](@ref rcAPI) |
| electromagnetic | TE and TM spherical modes | `SphericalModeTE`, `SphericalModeTM` | [spherical modes](@ref modesAPI) |
| electromagnetic | uniform static field | `UniformField` | [uniform static field](@ref uniformAPI) |
| acoustic | plane wave | `Acoustic.planeWave` | [acoustic plane wave](@ref ACpwAPI) |
| acoustic | monopole | `Acoustic.monopole` | [monopole](@ref ACpointAPI) |

The acoustic constructors are collected in the [`Acoustic` submodule](@ref ACsubmodule), so that their names do not clash with the electromagnetic ones.


---
## Defining a Scatterer

For all available scatterers a constructor with keyword arguments is available, e.g.,
```julia
sp = PECSphere(radius=1.0)                                        # electromagnetic
sp = DielectricSphere(radius=1.0, filling=Medium(4ε0, μ0))

sp = SoftSphere(radius=1.0)                                       # acoustic
sp = Spheroid{SoundHard}(equatorialRadius=1.0, polarRadius=0.5)   # an oblate spheroid
sp = Disc(SoundSoft; radius=1.0)
```
The physics of the scatterer has to match that of the excitation. An overview of all constructors and of the type hierarchy behind them is given in [Scatterers and Boundary Conditions](@ref scatterersConcept).


---
## Computing Fields

The incident, scattered, and total fields are computed by
```julia
F = field(ex, quantity)                # incident field: the excitation alone
F = scatteredfield(sp, ex, quantity)   # field scattered by the scatterer
F = field(sp, ex, quantity)            # total field: incident plus scattered
```
where the quantity wraps the locations, e.g.,
```julia
E  = field(sp, ex, ElectricField(point_cart))             # electromagnetic
H  = field(sp, ex, MagneticField(point_cart))
Φ  = field(sp, ex, ScalarPotential(point_cart))           # uniform static field
FF = scatteredfield(sp, ex, FarField(point_cart))

p  = field(sp, ex, Pressure(point_cart))                  # acoustic
γ₀ = field(sp, ex, PressureTrace(point_cart))             # on the surface
γ₁ = field(sp, ex, PressureNormalGradient(point_cart))    # on the surface
FF = scatteredfield(sp, ex, FarField(point_cart))
```
Which quantities are available for which excitation, and how the traces and far fields are defined, is described in [Quantities](@ref quantitiesConcept). All functions accept the keyword argument `parameter`, which controls the truncation of the series, see [Accuracy Settings](@ref).

!!! tip
    The locations are distributed over the available threads, for all excitations but the dipoles. Start Julia with several threads, e.g., `julia --threads=auto`, to make use of them.

!!! tip
    For a spheroid, the expensive modal coefficients can be computed once and reused for several quantities, see [Reusing the Modal Coefficients](@ref).


---
## Radar Cross Section

For an electromagnetic plane wave, the bistatic and the monostatic [radar cross section](@ref rcsPW) are computed by
```julia
σ = rcs(sp, ex, points_cart)   # bistatic

σ = rcs(sp, ex)                # monostatic
```


---
## Conversion Between Bases

Methods are provided to convert between Cartesian and spherical coordinates, whose convention is given under [Spherical Coordinates](@ref):

```julia
point_cart = SphericalScattering.sph2cart.(point_sph)

point_sph  = SphericalScattering.cart2sph.(point_cart)
```

Converting fields:

```julia
F_cart = SphericalScattering.convertSpherical2Cartesian.(F_sph,  point_sph)

F_sph  = SphericalScattering.convertCartesian2Spherical.(F_cart, point_sph)
```


---
## Plotting Fields

To visualize far-fields and near-fields the functions

```@docs
sphericalGridPoints
phiCutPoints
thetaCutPoints
```

as well as the functions

```julia
plotff(F, points_sph; scale="log", normalize=true, type="abs")

plotffcut(F, points; scale="log", normalize=true, format="polar")
```

are provided (after loading the [PlotlyJS](https://github.com/JuliaPlots/PlotlyJS.jl/tree/master) package). 
For more details see the [visualization of fields](@ref visualize) examples.
