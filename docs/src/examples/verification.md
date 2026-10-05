
# Code Verification

This package is well suited to verify the correctness of more involved numerical techniques to determine the scattering behavior of real-world objects, such as boundary element methods. Two examples are given below, both employing the packages [BEAST](https://github.com/krcools/BEAST.jl) and [CompScienceMeshes](https://github.com/krcools/CompScienceMeshes.jl): a PEC sphere and a sound-soft and a sound-hard sphere, each excited by a plane wave.


---
## Electromagnetic: PEC Sphere

The scattering from PEC objects can be determined by solving the electric field integral equation (EFIE) by the method of moments (MoM):

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
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide
# --- define points where the fields are computed
using SphericalScattering # use this package to get points on a spherical grid
points_cartFF, points_sphFF = sphericalGridPoints()
points_cartNF, points_sphNF = sphericalGridPoints(r=5.0)

# --- compute the fields 
EF_MoM = potential(MWSingleLayerField3D(𝑇), points_cartNF, u, RT)
HF_MoM = potential(BEAST.MWDoubleLayerField3D(wavenumber=κ), points_cartNF, u, RT) / (c * μ)
FF_MoM = -im * f / (2 * c) * potential(MWFarField3D(𝑇), points_cartFF, u, RT)
nothing # hide
end # hide
```

The fields for the same scenario can be obtained (more accurately) by this package:

```@example beast
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide
using SphericalScattering

sp = PECSphere(radius=spRadius)
ex = planeWave(frequency=f)

EF = scatteredfield(sp, ex, ElectricField(points_cartNF))
HF = scatteredfield(sp, ex, MagneticField(points_cartNF))
FF = scatteredfield(sp, ex, FarField(points_cartFF))
nothing #hide
end # hide
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

---
## Acoustic: Sound-Soft and Sound-Hard Sphere

The same approach applies to the acoustic scatterers. Both packages employ the time convention ``\mathrm{e}^{\mathrm{j}\omega t}``, so that the incident plane wave ``u^\mathrm{i} = \mathrm{e}^{-\mathrm{j} \kappa \hat{\bm e}_z \cdot \bm r}`` and the Green's function ``\mathrm{e}^{-\mathrm{j} \kappa R} / (4 \pi R)`` agree. As in the electromagnetic example, the sphere is embedded in the default medium [`Medium(ε0, μ0)`](@ref Medium), so that the wavenumber follows from the frequency as ``\kappa = 2 \pi f / c`` with ``c = 1 / \sqrt{\varepsilon_0 \mu_0}``:

```@example beastAcoustic
using BEAST, CompScienceMeshes
using SphericalScattering

# --- parameters
f = 7e7                 # frequency
c = 1 / sqrt(ε0 * μ0)   # wave speed in the default medium Medium(ε0, μ0)
κ = 2π * f / c          # wavenumber

# --- obtain a triangulation of the sphere
spRadius = 1.0 # radius of sphere
Γ = meshsphere(spRadius, 0.15)

# --- excitation by a plane wave travelling in z-direction
uⁱ = Helmholtz3D.planewave(; direction=ẑ, wavenumber=κ) # BEAST
ex = Acoustic.planeWave(; frequency=f)                  # this package

# --- points where the scattered fields are compared
points_cartNF, points_sphNF = sphericalGridPoints(r=5.0)
nothing # hide
```

On a sound-soft sphere the total pressure ``u = u^\mathrm{i} + u^\mathrm{s}`` vanishes. Its Neumann trace ``\sigma = \partial_n u`` solves the first-kind equation ``\mathcal{S} \sigma = u^\mathrm{i}`` with the single-layer operator ``\mathcal{S}``, and it radiates the scattered field ``u^\mathrm{s} = -\mathcal{S}[\sigma]``:

```@example beastAcoustic
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide
# --- sound-soft: piecewise constant Neumann trace
X = lagrangecxd0(Γ)
S = assemble(Helmholtz3D.singlelayer(; wavenumber=κ), X, X)
σ_MoM = S \ assemble(strace(uⁱ, Γ), X)
us_soft_MoM = -potential(HH3DSingleLayerNear(; wavenumber=κ), points_cartNF, σ_MoM, X; type=ComplexF64)

# --- the same quantities by this package, the trace at the centroids of the triangles
sp = SoftSphere(radius=spRadius)
σ = field(sp, ex, PressureNormalGradient(X.pos))
us_soft = scatteredfield(sp, ex, Pressure(points_cartNF))
nothing # hide
end # hide
```

On a sound-hard sphere the normal derivative of the total pressure vanishes. Its Dirichlet trace ``u`` solves the first-kind equation ``\mathcal{W} u = \partial_n u^\mathrm{i}`` with the hypersingular operator ``\mathcal{W}``, and it radiates the scattered field ``u^\mathrm{s} = \mathcal{D}[u]`` by the double-layer potential:

```@example beastAcoustic
begin # hide
redirect_stderr(devnull) # the progress bars would clutter the output # hide
# --- sound-hard: piecewise linear Dirichlet trace
Y = lagrangec0d1(Γ)
W = assemble(Helmholtz3D.hypersingular(; wavenumber=κ), Y, Y)
u_MoM = W \ assemble(∂n(uⁱ), Y)
us_hard_MoM = potential(HH3DDoubleLayerNear(; wavenumber=κ), points_cartNF, u_MoM, Y; type=ComplexF64)

# --- the same quantities by this package, the trace at the vertices
sp = HardSphere(radius=spRadius)
u = field(sp, ex, PressureTrace(Y.pos))
us_hard = scatteredfield(sp, ex, Pressure(points_cartNF))
nothing # hide
end # hide
```

The centroids and the vertices of the mesh do not lie exactly on the sphere; as the traces depend only on the direction of a location, they are evaluated on the sphere nevertheless, see [Quantities](@ref quantitiesConcept). The agreement is again determined as a worst case relative error of all evaluated points:

```@example beastAcoustic
# --- relative worst case errors in percent
relerr(a, b) = round(maximum(abs.(a - b)) / maximum(abs.(b)) * 100, digits=2)

print("sound-soft, Neumann trace error:   $(relerr(σ_MoM, σ)) %\n")
print("sound-soft, scattered field error: $(relerr(us_soft_MoM, us_soft)) %\n")
print("sound-hard, Dirichlet trace error: $(relerr(u_MoM, u)) %\n")
print("sound-hard, scattered field error: $(relerr(us_hard_MoM, us_hard)) %\n")
```

!!! note
    The first-kind equations are not uniquely solvable at the interior resonances of the sphere: the single-layer operator fails at the wavenumbers of the interior Dirichlet problem, the first of which is ``\kappa a = \pi``, the hypersingular operator at those of the interior Neumann problem, the first of which is ``\kappa a \approx 2.08``. The frequency above is chosen such that ``\kappa a \approx 1.47`` lies below both.


!!! tip
    When [testing this package](@ref tests), the packages [BEAST](https://github.com/krcools/BEAST.jl) and 
    [CompScienceMeshes](https://github.com/krcools/CompScienceMeshes.jl) are used to define several functional tests.
