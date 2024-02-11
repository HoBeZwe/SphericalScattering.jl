
# SphericalScattering.jl

This package provides semi-analytical solutions to the scattering of time harmonic and static electromagnetic fields from spherical objects (amongst others known as Mie solutions or Mie scattering). 
To this end, series expansions are evaluated. Special care is taken to obtain accurate solutions down to the static limit.

!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` and SI units are used everywhere.


---
## Installation

Installing SphericalScattering is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add SphericalScattering 
```


---
## Overview

The following aspects are implemented (✔) and planned (⌛):

##### Available incident fields:
- ✔ Plane wave
- ✔ Field of electric/magnetic ring current
- ✔ Field of electric/magnetic dipole
- ✔ TE/TM spherical vector waves
- ✔ Uniform static electric field
- ⌛ Static charge(s)

##### Available scattering objects:
- ✔ PEC sphere
- ⌛ PMC sphere
- ⌛ Dielectric sphere 
- ⌛ Multilayer dielectric sphere 
- ⌛ Multilayer dielectric sphere with PEC core 
- ✔ Dielectric sphere with thin impedance layer

##### Available quantities (where applicable):
- ✔ Far-fields
- ✔ Near-fields (electric & magnetic)
- ✔ Radar cross section (RCS)
- ⌛ Surface currents
- ✔ Scalar potentials 
- ✔ Displacement fields 
- ✔ Scalar potential jump 

---
##### Detailed implementation status:

| spheres                              | plane wave | el. ring current | mag. ring current | el. dipole | mag. dipole | TE/TM modes | uniform static field | static charge(s) |
|--------------------------------------|------------|------------------|-------------------|------------|-------------|-------------|----------------------|------------------|
| PEC                                  |      ✔     |        ✔         |         ✔         |      ✔     |       ✔     |      ✔      |           ✔          |        ⌛         |
| PMC                                  |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ⌛          |        ⌛        |
| Dielectric                           |      ✔     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Multilayer dielectric                |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Multilayer dielectric with PEC core  |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Dielectric with thin impedance layer |      ➖     |        ➖         |         ➖         |      ➖     |       ➖    |      ➖      |           ✔          |        ➖        |




---
## [References](@id refs)

The implementation is based on
- [1] J.-M. Jin, *Theory and Computation of Electromagnetic Fields*, Second edition. Hoboken, New Jersey: John Wiley & Sons, Inc, 2015.
- [2] G. T. Ruck, D. E. Barrick, W. D. Stuart, C. K. Krichbaum, *Radar Cross Section Handbook*, Volume 1, New York: Plenum Press, 1970.
- [3] Sihvola, Ari & Lindell, Ismo., *Transmission Line Analogy for Calculating the Effective Permittivity of Mixtures with Spherical Multilayer Scatterers*, Journal of Electromagnetic Waves and Applications, Volume 2, Pages 741-756, 1988.
- [4] J. D. Jackson, *Classical Electrodynamics*, New York: Wiley, 3rd ed., 1999.
- [5] J. E. Hansen, ed., *Spherical Near-field Antenna Measurements*, The Institution of Engineering and Technology, Michael Faraday House, Six Hills Way, Stevenage SG1 2AY, UK: IET, 1988.
- [6] T. B. Jones, Ed., *Models for layered spherical particles*, in Electromechanics of Particles, Cambridge: Cambridge University Press, 1995, pp. 227–235. doi: 10.1017/CBO9780511574498.012.

!!! note
    If you use this software, please cite our [JOSS article](https://doi.org/10.21105/joss.05820):

    B. Hofmann, P. Respondek, and S. B. Adrian, *Sphericalscattering: a Julia package for electromagnetic scattering from spherical objects*, Journal of Open Source Software, vol. 8, no. 91, Nov. 2023, doi: 10.21105/joss.05820.