
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

!!! tip
    When [testing this package](@ref tests), the packages [BEAST](https://github.com/krcools/BEAST.jl) and 
    [CompScienceMeshes](https://github.com/krcools/CompScienceMeshes.jl) are used to define several functional tests.


---
## [Radar Cross Section](@id rcsApplication) 

The monostatic and the bistatic [radar cross section (RCS)](@ref rcsPW) can be evaluated.

### Monostatic RCS

The monostatic radar cross section as a function of the sphere radius can be computed as follows (compare also the plot in [jinTheoryComputationElectromagnetic2015; pp. 352ff](@cite)):
```@example RCS
using SphericalScattering
using PlotlyJS

f = 1e8
c = 2.99792458e8    # speed of light
λ = c / f

RCS = Float64[]
aTλ = Float64[]

# --- compute RCS
for rg in 0.01:0.01:3.0

    a = λ*rg

    monoRCS = rcs(PECSphere(radius=a), planeWave(frequency=f))

    push!(RCS, monoRCS / (π * a^2))
    push!(aTλ, rg)
end

# --- plot
layout = Layout(
    yaxis=attr(title_text="RCS / πa² in dB"),
    xaxis=attr(title_text="a / λ")
)

plot(scatter(; x=aTλ, y=10*log10.(RCS), mode="lines+markers"), layout)
t = plot(scatter(; x=aTλ, y=10*log10.(RCS), mode="lines+markers"), layout) # hide
savefig(t, "monoRCS.html"); nothing # hide
```

```@raw html
<object data="../monoRCS.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

### Bistatic RCS

The bistatic radar cross section along a ϑ-cut can be computed as follows (compare also the plot in [jinTheoryComputationElectromagnetic2015; pp. 351ff](@cite)):
```@example RCS
using StaticArrays

# --- points
ϑ = [i*π/500 for i in 0:500]
φ = 0
points_cart = [SphericalScattering.sph2cart(SVector(1.0, ϑi, φ)) for ϑi in ϑ]

# --- compute RCS
biRCS = rcs(PECSphere(radius=λ), planeWave(frequency=1e8), points_cart) / λ^2

# --- plot
layout = Layout(
    yaxis=attr(title_text="RCS / λ² in dB"),
    xaxis=attr(title_text="ϑ in degree")
)

plot(scatter(; x=ϑ*180/π, y=10*log10.(biRCS), mode="lines+markers"), layout)
t = plot(scatter(; x=ϑ*180/π, y=10*log10.(biRCS), mode="lines+markers"), layout) # hide
savefig(t, "biRCS.html"); nothing # hide
```

```@raw html
<object data="../biRCS.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```



---
## [Visualization of Fields](@id visualize) 

This package provides several means to directly visualize quantities of the scattering setup.

!!! note
    In order to make the plotting functionality available the package [PlotlyJS](https://github.com/JuliaPlots/PlotlyJS.jl/tree/master) has to be loaded.
    (It is a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).)


### Plotting Far-Field Patterns

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
<object data="../plotPatternHDPEC.html" type="text/html"  style="width:100%;height:50vh;"> </object>
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
<object data="../plotPatternHDtot.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```



### Plotting Far-Field Cuts

Sphercial cuts can be conveniently obtained by:

```@example plotcuts
using SphericalScattering
using LinearAlgebra, StaticArrays

# --- excite PEC sphere by magnetic ring current
orient = normalize(SVector(0.0,1.0,1.0))
ex = magneticRingCurrent(frequency=1e8, orientation=orient, center=2*orient, radius=0.2)

sp = PECSphere(radius = 1.0)

# --- evaluate fields at φ = 5° cut
points_cart, points_sph = phiCutPoints(5) # analogously, thetaCutPoints can be used

FF = scatteredfield(sp, ex, FarField(points_cart))
nothing # hide
```

The cut can then be plotted as:

```@example plotcuts
using PlotlyJS
plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="polar")
t = plotffcut(norm.(FF), points_sph, normalize=true, scale="log", format="polar") # hide
savefig(t, "plotcut.html"); nothing # hide
```

```@raw html
<object data="../plotcut.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```


### Plotting Near-Field Cuts

The near-field of an electric ring current can, e.g., be visualized in the xz-plane as:

```@example heatmap
using SphericalScattering
using LinearAlgebra, StaticArrays

f = 1e8             # frequency
c = 2.99792458e8    # speed of light
λ = c / f           # wavelength

ex = electricRingCurrent(frequency=1e8, center=SVector(0.,0,0), radius=3*λ)

# --- define points in the xz plane
res = λ/15

points_cart = [SVector(x, 0.0, z) for z in -5λ:res:5λ, x in -5λ:res:5λ]
points_sph = SphericalScattering.cart2sph.(points_cart)

# --- evaluate the fields
E = field(ex, ElectricField(points_cart))

nothing # hide
```

To plot the fields a heatmap can be employed:

```@example heatmap
using PlotlyJS

# --- convert data to spherical coordinates
Esph = SphericalScattering.convertCartesian2Spherical.(E, points_sph)
Eφ = [Esph[i][3] for (i, j) in enumerate(Esph)] # extract φ-component

# --- truncate large values (for plot)
data = abs.(Eφ)
data[data .> 250] .= 250

# --- plot heatmap
layout = Layout(
    yaxis=attr(title_text="y/λ"),
    xaxis=attr(title_text="x/λ")
)

plot(heatmap(x=-5:1/15:5,y=-5:1/15:5,z=data, colorscale="Jet"), layout)
t = plot(heatmap(x=-5:1/15:5,y=-5:1/15:5,z=data, colorscale="Jet"), layout) # hide
savefig(t, "plotNF.html"); nothing # hide
```

```@raw html
<object data="../plotNF.html" type="text/html"  style="width:60%;height:50vh;"> </object>
```

Or instead of the magnitude a snapshot can be plotted

```@example heatmap
data = real.(Eφ)
data[data .> 250] .= 250

# --- plot heatmap
layout = Layout(
    yaxis=attr(title_text="y/λ"),
    xaxis=attr(title_text="x/λ")
)

plot(heatmap(x=-5:1/15:5,y=-5:1/15:5,z=data, colorscale="Jet"), layout)
t = plot(heatmap(x=-5:1/15:5,y=-5:1/15:5,z=data, colorscale="Jet"), layout) # hide
savefig(t, "plotNF2.html"); nothing # hide
```

```@raw html
<object data="../plotNF2.html" type="text/html"  style="width:60%;height:50vh;"> </object>
```