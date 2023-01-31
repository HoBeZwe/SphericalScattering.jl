
# SphericalScattering.jl

This package provides semi-analytical solutions to the scattering of time harmonic and static electromagnetic fields from spherical objects. 
To this end, series expansions are evaluated. Special care is taken to obtain accurate solutions down to the static limit.


---
## Overview

The following scenarios are implemented (✓) and planned (⌛):

##### Available incident fields:
- ✓ Plane wave
- ✓ Field of electric/magnetic ring current
- ✓ Field of electric/magnetic dipole
- ✓ TE/TM spherical vector waves
- ✓ Uniform static electric field

##### Available scattering objects:
- ✓ PEC sphere
- ⌛ PMC sphere
- ⌛ Dielectric sphere (so far only for uniform field and the plane-wave)
- ⌛ Multilayer dielectric sphere (so far only for uniform field)
- ⌛ Multilayer dielectric sphere with PEC core (so far only for uniform field)


!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` is used everywhere.

---
## Installation

Installing SphericalScattering is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add SphericalScattering 
```

---
## [References](@id refs)

The implementation is based on
- [1] J.-M. Jin, *Theory and Computation of Electromagnetic Fields*, Second edition. Hoboken, New Jersey: John Wiley & Sons, Inc, 2015.
- [2] G. T. Ruck, D. E. Barrick, W. D. Stuart, C. K. Krichbaum, *Radar Cross Section Handbook*, Volume 1, New York: Plenum Press, 1970.
- [3] Sihvola, Ari & Lindell, Ismo., *Transmission Line Analogy for Calculating the Effective Permittivity of Mixtures with Spherical Multilayer Scatterers*, Journal of Electromagnetic Waves and Applications, Volume 2, Pages 741-756, 1988.
- [4] J. D. Jackson, *Classical Electrodynamics*, New York: Wiley, 3rd ed., 1999.
- [5] J. E. Hansen, ed., *Spherical Near-field Antenna Measurements*, The Institution of Engineering and Technology, Michael Faraday House, Six Hills Way, Stevenage SG1 2AY, UK: IET, 1988.
- [6] T. B. Jones, Ed., *Models for layered spherical particles*, in Electromechanics of Particles, Cambridge: Cambridge University Press, 1995, pp. 227–235. doi: 10.1017/CBO9780511574498.012.
