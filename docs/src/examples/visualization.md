
# [Visualization of Fields](@id visualize)

This package provides several means to directly visualize quantities of the scattering setup.

!!! note
    In order to make the plotting functionality available the package [PlotlyJS](https://github.com/JuliaPlots/PlotlyJS.jl/tree/master) has to be loaded.
    (It is a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).)


## Plotting Far-Field Patterns

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
<object data="../../examples/plotPatternHDPEC.html" type="text/html"  style="width:100%;height:50vh;"> </object>
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
<object data="../../examples/plotPatternHDtot.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```



## Plotting Far-Field Cuts

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
<object data="../../examples/plotcut.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```


## Plotting Near-Field Cuts

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
<object data="../../examples/plotNF.html" type="text/html"  style="width:60%;height:50vh;"> </object>
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
<object data="../../examples/plotNF2.html" type="text/html"  style="width:60%;height:50vh;"> </object>
```
