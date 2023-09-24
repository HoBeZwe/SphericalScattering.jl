
# Application Examples

Two classical application examples are discussed below.

---
## Code Verification

This package is well suited to verify the correctness of more involved numerical techniques to determine the scattering behavior of real-world objects.
For example, the scattering from PEC objects can be determined by solving the electric field integral equation (EFIE) by the method of moments (MoM).
To this end, the packages [BEAST](https://github.com/krcools/BEAST.jl) and [CompScienceMeshes](https://github.com/krcools/CompScienceMeshes.jl) can be employed:

```@example beast
using BEAST, CompScienceMeshes

# --- parameters
f = 1e8             # frequency
c = 2.99792458e8    # speed of light
μ = 4π * 1e-7       # permeability
κ = 2π * f / c      # wavenumber

# --- obtain a triangulation of the sphere
spRadius = 1.0 # radius of sphere
Γ = meshsphere(spRadius, 0.25)

# --- define basis functions on the triangulation
RT = raviartthomas(Γ)

# --- excitation by plane wave
𝐸 = Maxwell3D.planewave(; direction=ẑ, polarization=x̂, wavenumber=κ)
𝑒 = n × 𝐸 × n

# --- integral operator 
𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)

# --- assemble matrix and RHS of the LSE
e = -assemble(𝑒, RT)
T = assemble(𝑇, RT, RT)

# --- solve the LSE 
u = T \ e
```

Now that we know the expansion coefficients for the basis functions we can compute the scattered fields:


```@example beast
# --- define points where the fields are computed
using SphericalScattering # use this package to get points on a spherical grid
points_cartFF, points_sphFF = sphericalGridPoints()
points_cartNF, points_sphNF = sphericalGridPoints(r=5.0)

# --- compute the fields 
EF_MoM = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
HF_MoM = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), points_cartNF, u, RT) / (c * μ)
FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)
nothing # hide
```

The fields for the same scenario can be obtained (more accurately) by this package:

```@example beast
using SphericalScattering

sp = PECSphere(radius=spRadius)
ex = planeWave(frequency=f)

EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
FF = scatteredfield(sp, ex, FarField(points_cartFF))
nothing #hide
```

The agreement between the two solutions can be determined as a worst case relative error of all evaluated points:

```@example beast
using LinearAlgebra

# --- relative worst case errors in percent
diff_EF = round(maximum(norm.(EF - EF_MoM) ./ maximum(norm.(EF))) * 100, digits=2)
diff_HF = round(maximum(norm.(HF - HF_MoM) ./ maximum(norm.(HF))) * 100, digits=2)
diff_FF = round(maximum(norm.(FF - FF_MoM) ./ maximum(norm.(FF))) * 100, digits=2)

print("E-field error: $diff_EF %\n")
print("H-field error: $diff_HF %\n")
print("Far-field error: $diff_FF %\n")
```

!!! note
    When [testing this package](@ref tests), the packages [BEAST](https://github.com/krcools/BEAST.jl) and 
    [CompScienceMeshes](https://github.com/krcools/CompScienceMeshes.jl) are used to define several functional tests.




---
## Visualization of Scattering 

This package provides several means to directly visualize quantities of the scattering setup.

!!! note
    In order to make the plotting functionality available the package [PlotlyJS](https://github.com/JuliaPlots/PlotlyJS.jl/tree/master) has to be loaded.
    (It is a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).)


### Plot the Far-Field Patterns

As an example consider a Hertzian dipole that excites a PEC sphere:

```@example plotpattern
using SphericalScattering
using LinearAlgebra, StaticArrays

# --- excite PEC sphere by Hertzian dipole
orient = normalize(SVector(0.0,1.0,1.0))
ex = HertzianDipole(frequency=1e8, orientation=orient, position=2*orient)

sp = PECSphere(radius = 1.0)

# --- evaluate fields at spherical grid points
points_cart, points_sph = sphericalGridPoints()

FF = scatteredfield(sp, ex, FarField(points_cart))
nothing # hide
```

The 3D pattern can then be evaluated as

```@example plotpattern
using PlotlyJS
plotff(FF, points_sph, scale="linear", normalize=true, type="abs")
t = plotff(FF, points_sph, scale="linear", normalize=true, type="abs") # hide
savefig(t, "plotPatternHDPEC.html"); nothing # hide
```

```@raw html
<object data="plotPatternHDPEC.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

Alternatively, the field of the dipole itself or the total field can be plotted:

```@example plotpattern
FF = field(sp, ex, FarField(points_cart)) # total field

# --- plot only the φ-component in logarithmic scale   
plotff(FF, points_sph, scale="log", normalize=true, type="phi") 
t = plotff(FF, points_sph, scale="log", normalize=true, type="phi") # hide
savefig(t, "plotPatternHDtot.html"); nothing # hide
```

```@raw html
<object data="plotPatternHDtot.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```