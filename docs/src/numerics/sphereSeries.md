
# [Series for Spheres](@id sphereSeries)

The fields of the spheres are evaluated as series in the spherical wave functions, the Mie series [jinTheoryComputationElectromagnetic2015, bowmanElectromagneticAcousticScattering1970](@cite). Each term of degree ``n`` combines a coefficient of the incident field, a scattering coefficient, which depends on the sphere but not on the excitation, a spherical Bessel or Hankel function of the radial coordinate, and Legendre functions of the angles. The series of the electromagnetic and of the acoustic spheres are discussed below, as well as how they are evaluated.

Two excitations require no series: the field of a [uniform static field](@ref uniformEx) is given in closed form, and a [spherical mode](@ref modesAPI) excites a single term.


---
## [Symmetry](@id sphereSymmetry)

Each series is formulated in a frame in which the excitation is as symmetric as possible:

- **Electromagnetic.** The plane wave travels along ``\hat{\bm e}_z`` with polarization ``\hat{\bm e}_x``; the dipoles and ring currents are oriented along the ``z``-axis. An excitation of arbitrary orientation is related to this frame by a rotation of the locations and of the fields, see [Rotations](@ref rotationDetails).
- **Acoustic.** The scattered pressure depends only on the angle between the observation direction and the symmetry axis of the excitation, which is the direction of incidence of a plane wave and the line from the center of the sphere to a monopole. The series is a sum over Legendre polynomials of the cosine of this angle, so that no rotation is needed.


---
## [Electromagnetic Spheres](@id emSphereSeries)

The series of the PEC and the dielectric sphere follow [jinTheoryComputationElectromagnetic2015](@cite): pp. 347ff for the plane wave, pp. 368ff for the ring currents, and a generalization of the analysis on pp. 374ff for the dipoles. The magnetic dipole and ring current follow from the electric ones by the [duality relations](@ref dualityRelations). Each series is formulated in the frame described under [Symmetry](@ref sphereSymmetry).

#### Spherical Modes

A spherical mode, see its [definition](@ref modesDefinition), is scattered by a PEC sphere as a single mode. Matching incoming and outgoing waves to fulfill the boundary condition ``\bm e_\mathrm{tan} = \bm 0`` yields the scattering coefficients
```math
\xi_\mathrm{TE} = -\cfrac{\mathrm{H}^{(2)}_{n + 0.5}(k a)}{\mathrm{H}^{(1)}_{n + 0.5}(k a)}
```
and
```math
\xi_\mathrm{TM} = -\cfrac{\mathrm{H}'^{(2)}_{n + 0.5}(k a)}{\mathrm{H}'^{(1)}_{n + 0.5}(k a)}
```
where ``\mathrm{H}^{(\nu)}_{n}(x)`` denotes the Hankel function of ``\nu``-th kind and ``n``-th order.

The scattered fields ``\bm e^\mathrm{sc}`` are then given by
```math
\bm e_\mathrm{TE}^\mathrm{sc} = \xi_\mathrm{TE} k \sqrt{Z_\mathrm{F}} \bm{f}_{1mn}^{(2)}
```
and
```math
\bm e_\mathrm{TM}^\mathrm{sc} = \xi_\mathrm{TM} k \sqrt{Z_\mathrm{F}} \bm{f}_{2mn}^{(2)} \,.
```


---
## [Acoustic Spheres](@id acSphereSeries)

The series of the [sound-hard and sound-soft spheres](@ref acScattererAPI) follow [bowmanElectromagneticAcousticScattering1970](@cite). In contrast to the electromagnetic case they start at ``n=0``: the monopole term contributes and dominates the low-frequency limit.

#### Plane Wave

Expanding the incident pressure in the spherical waves about the center of the sphere,
```math
p_\mathrm{PW}(\bm r) = a \sum_{n=0}^\infty (2n+1) (-\mathrm{j})^n j_n(k r) P_n(\cos \vartheta) \,,
```
where ``\vartheta`` is measured from the direction of incidence ``\hat{\bm d}`` and ``P_n`` denotes the Legendre polynomials, the scattered pressure follows as
```math
p^\mathrm{sc}(\bm r) = a \sum_{n=0}^\infty (2n+1) (-\mathrm{j})^n b_n h_n^{(2)}(k r) P_n(\cos \vartheta) \,.
```
Applying the boundary condition term by term yields the scattering coefficients
```math
b_n = -\cfrac{j_n^\prime(k r_\mathrm{s})}{h_n^{(2)\prime}(k r_\mathrm{s})}
\qquad \text{and} \qquad
b_n = -\cfrac{j_n(k r_\mathrm{s})}{h_n^{(2)}(k r_\mathrm{s})}
```
for a sound-hard and a sound-soft sphere of radius ``r_\mathrm{s}``, respectively: on a sound-hard surface the normal velocity, and hence the radial derivative of the total pressure, vanishes, whereas on a sound-soft (pressure release) surface the total pressure vanishes [bowmanElectromagneticAcousticScattering1970](@cite).

The far field is defined as
```math
p^\mathrm{sc}_\infty(\hat{\bm r}) = \lim_{r \rightarrow \infty} r \, \mathrm{e}^{\mathrm{j} k r} p^\mathrm{sc}(\bm r)
                                  = \cfrac{\mathrm{j} a}{k} \sum_{n=0}^\infty (2n+1) b_n P_n(\cos \vartheta) \,,
```
that is, the factor ``\mathrm{e}^{-\mathrm{j} k r} / r`` is omitted, as for the electromagnetic excitations. Note that the order-dependent factors cancel, since ``(-\mathrm{j})^n \mathrm{j}^{n+1} = \mathrm{j}`` holds for every ``n``.

#### Monopole

The addition theorem of the free-space Green's function,
```math
\cfrac{\mathrm{e}^{-\mathrm{j} k |\bm r - \bm r_0|}}{4 \pi |\bm r - \bm r_0|}
    = \cfrac{-\mathrm{j} k}{4 \pi} \sum_{n=0}^\infty (2n+1) j_n(k r_<) h_n^{(2)}(k r_>) P_n(\cos \vartheta) \,,
```
with ``r_< = \min(r, r_0)``, ``r_> = \max(r, r_0)``, ``r_0 = |\bm r_0|``, the Legendre polynomials ``P_n`` and ``\vartheta`` measured from the direction ``\hat{\bm r}_0`` towards the monopole, expands the incident pressure in the spherical waves about the center of the sphere. Since the monopole lies outside the sphere, ``r_> = r_0`` holds on its surface, where the boundary condition is imposed. The scattered pressure therefore reads
```math
p^\mathrm{sc}(\bm r) = \cfrac{-\mathrm{j} k a}{4 \pi} \sum_{n=0}^\infty (2n+1) h_n^{(2)}(k r_0) \, b_n \, h_n^{(2)}(k r) P_n(\cos \vartheta)
```
with the very same scattering coefficients ``b_n`` as for the plane wave. Applying the boundary condition term by term leaves the radial structure of the series untouched, hence the scattering coefficients do not depend on the excitation; only the coefficients of the incident expansion do.

!!! tip
    The series converges like ``(r_\mathrm{s}^2 / (r_0 r))^n``. A monopole very close to the surface therefore requires many terms, which may exceed the numerical range of the spherical Hankel functions before the series converges; a message is printed should this happen.

#### Neumann Trace

For a normal ``\hat{\bm n} \neq \hat{\bm r}`` the Neumann trace of the scattered field picks up the tangential part of the gradient as well:
```math
\hat{\bm n} \cdot \nabla p^\mathrm{sc} = (\hat{\bm n} \cdot \hat{\bm r}) \cfrac{\partial p^\mathrm{sc}}{\partial r}
    + (\hat{\bm n} \cdot \hat{\bm \vartheta}) \cfrac{1}{r} \cfrac{\partial p^\mathrm{sc}}{\partial \vartheta} \,.
```
Both contributions are included, so that providing normals is supported here as well.


---
## Special Functions

The spherical Bessel and Hankel functions are obtained from the cylindrical ones of half-integer order,
```math
j_n(x) = \sqrt{\cfrac{\pi}{2x}} \, J_{n + 1/2}(x) \,, \qquad
h^{(2)}_n(x) = \sqrt{\cfrac{\pi}{2x}} \, H^{(2)}_{n + 1/2}(x) \,,
```
which are provided by [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl). Their derivatives, and those of the Riccati-Bessel functions ``x j_n(x)`` and ``x h^{(2)}_n(x)``, follow from recurrence relations in the order. The Legendre functions are computed by their three-term recurrence in the degree, or by [LegendrePolynomials.jl](https://github.com/jishnub/LegendrePolynomials.jl).


---
## Truncation

In contrast to the [spheroidal series](@ref ACspheroidAccuracy), whose truncation is fixed before any field is evaluated, the series of the spheres are summed at each location until they have converged: terms are added until the relative contribution of the last one drops below the `relativeAccuracy` of [`Parameter`](@ref), which defaults to ``10^{-12}``, but at least ten terms are taken. The number of terms thus adapts by itself to the size of the sphere and to the distance of the location.

!!! note
    The `nmax` of [`Parameter`](@ref) has no effect on the series of the spheres; it fixes the degree of the spheroidal series only.

A term of one parity may vanish identically, e.g., a Legendre polynomial of odd degree in the plane perpendicular to the symmetry axis. Such a term must not be mistaken for convergence: depending on the series, the relative contribution is therefore measured only for the terms which do not vanish, or only for every second degree.


---
## Large Degrees and the Static Limit

For a small argument the spherical Hankel functions grow rapidly with the degree, roughly like ``(2n-1)!! / x^{n+1}``, while the scattering coefficients decay accordingly, so that their products stay small. At low frequencies, or for a large number of terms, the Hankel functions alone may exceed the range of the floating-point numbers before the series has converged. The summation then stops with the terms computed so far: the acoustic series print a message (`did not converge: n=…`) when this happens, the electromagnetic ones stop silently.

In the cases checked, the series remain accurate down to the static limit. For instance, at ``ka \approx 6 \cdot 10^{-8}`` the scattered pressure of a sound-soft sphere equals ``-a / r`` times the incident one to all printed digits, and the scattered electric field of a PEC sphere at a frequency of ``1\,\mathrm{mHz}`` is that of the induced static dipole.


---
## Validation

The series are checked by the [tests](@ref tests) of the package: the time-harmonic electromagnetic solutions against boundary element solutions computed by [BEAST](https://github.com/krcools/BEAST.jl), and the acoustic ones by their boundary conditions, the reciprocity of source and observation point, the far-field limit, and their limiting cases, such as a distant monopole approaching a plane wave.
