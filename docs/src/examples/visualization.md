
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



## [Plotting Surface Densities](@id visualizeSurface)

The densities on the surface of an acoustic scatterer can be plotted on the scatterer itself, e.g., the total pressure on a tilted, sound-hard oblate spheroid in air, excited by an oblique plane wave:

```@example surface
using SphericalScattering, StaticArrays, LinearAlgebra
using PlotlyJS
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide

air = Medium(1 / (1.2 * 343.0^2), 1.2)
ex = Acoustic.planeWave(frequency=500.0, embedding=air, direction=normalize(SVector(1.0, 0.0, 1.0)))

sp = Spheroid{SoundHard}(equatorialRadius=0.5, polarRadius=0.2, axis=normalize(SVector(0.0, 0.3, 1.0)))

p = plotsurface(sp, ex, PressureTrace)
savefig(p, "plotSurface.html"); rm("fort.40"; force=true); nothing # hide
end # hide
```

```@raw html
<object data="../../examples/plotSurface.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

The two faces of a disc are plotted apart. For a sound-soft disc, the Neumann trace is singular at the rim, so that the colour scale is clipped at a quantile of the magnitudes:

```@example surface
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide
d = Disc(SoundSoft; radius=0.5)

p = plotsurface(d, ex, PressureNormalGradient; type="real", clip=0.95, gap=0.1)
savefig(p, "plotSurfaceDisc.html"); rm("fort.40"; force=true); nothing # hide
end # hide
```

```@raw html
<object data="../../examples/plotSurfaceDisc.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```
