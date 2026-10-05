
# [Radar Cross Section](@id rcsPW)

To compute the bistatic radar cross section (RCS) [jinTheoryComputationElectromagnetic2015; pp. 350ff](@cite)
```math
\sigma (\vartheta, \varphi) = \lim_{r\rightarrow \infty} \left( 4 \pi r^2 \frac{{|e^\mathrm{sc}|}^2}{{|e^\mathrm{inc}|}^2} \right)
```
the function
```julia
σ = rcs(sp, ex, points_cart)
```
is provided. For the monostatic RCS, the function
```julia
σ = rcs(sp, ex)
```
is provided.

!!! note
    The RCS is (so far) only defined for a plane wave excitation.


---
## [Examples](@id rcsApplication) 

The monostatic and the bistatic radar cross section can be evaluated.

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
<object data="../../electromagnetic/monoRCS.html" type="text/html"  style="width:100%;height:50vh;"> </object>
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
<object data="../../electromagnetic/biRCS.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```
