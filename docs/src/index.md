
# SphericalScattering.jl

This package provides semi-analytical solutions to the scattering of time harmonic electromagnetic fields from spherical objects. 
To this end, series expansions are evaluated. Special care is taken to obtain accurate solutions down to the static limit.


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
- ✓/⌛ Dielectric sphere (only for uniform field)
- ✓/⌛ Mulitlayer dielectric sphere (only for uniform field)
- ✓/⌛ Mulitlayer dielectric sphere with PEC core (only for uniform field)


!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` is used.

## Installation

Installing SphericalScattering is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add SphericalScattering 
```

## References

The implementation is based on
- J.-M. Jin, Theory and Computation of Electromagnetic Fields, Second edition. Hoboken, New Jersey: John Wiley & Sons, Inc, 2015.
- G. T. Ruck, D. E. Barrick, W. D. Stuart, C. K. Krichbaum, Radar Cross Section Handbook, Volume 1, New York: Plenum Press, 1970.
- Sihvola, Ari & Lindell, Ismo. (1988). Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers. Journal of Electromagnetic Waves and Applications - J ELECTROMAGNET WAVE APPLICAT. 2. 741-756. 