
# SphericalScattering.jl

This package provides semi-analytical solutions to the scattering of time-harmonic and static electromagnetic fields as well as time-harmonic acoustic fields from spherical objects (amongst others known as Mie solutions or Mie scattering). 
To this end, series expansions are evaluated. Special care is taken to obtain accurate solutions down to the static limit.

!!! note
    A time convention of ``\mathrm{e}^{\,\mathrm{j}\omega t}`` and SI units are used everywhere.

!!! note
    If you use this software, please cite our [JOSS article](https://doi.org/10.21105/joss.05820):

    B. Hofmann, P. Respondek, and S. B. Adrian, *Sphericalscattering: a Julia package for electromagnetic scattering from spherical objects*, Journal of Open Source Software, vol. 8, no. 91, Nov. 2023, doi: 10.21105/joss.05820.


---
## Installation

Installing SphericalScattering is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add SphericalScattering 
```


---
## Feature Overview

The following aspects are implemented (✔) and planned (⌛):

- Electromagnetic

| spheres                              | plane wave | el. ring current | mag. ring current | el. dipole | mag. dipole | TE/TM modes | uniform static field | static charge(s) |
|--------------------------------------|------------|------------------|-------------------|------------|-------------|-------------|----------------------|------------------|
| PEC                                  |      ✔     |        ✔         |         ✔         |      ✔     |       ✔     |      ✔      |           ✔          |        ⌛         |
| PMC                                  |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ⌛          |        ⌛        |
| Dielectric                           |      ✔     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Multilayer dielectric                |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Multilayer dielectric with PEC core  |      ⌛     |        ⌛         |         ⌛         |      ⌛     |       ⌛     |      ⌛      |           ✔          |        ⌛        |
| Dielectric with thin impedance layer |      ➖     |        ➖         |         ➖         |      ➖     |       ➖    |      ➖      |           ✔          |        ➖        |

- Acoustic

| objects                              | plane wave | monopole |
|--------------------------------------|------------|----------|
| Sphere sound-hard                    |      ✔     |     ✔    | 
| Sphere sound-soft                    |      ✔     |     ✔    | 
| Prolate Spheroid sound-hard          |      ⌛     |     ⌛    |
| Prolate Spheroid sound-soft          |      ⌛     |     ⌛    |
| Oblate Spheroid sound-hard           |      ⌛     |     ⌛    |
| Oblate Spheroid sound-soft           |      ⌛     |     ⌛    |
| Disc sound-hard                      |      ⌛     |     ⌛    |
| Disc sound-soft                      |      ⌛     |     ⌛    | 


---
##### Available incident fields:

- Electromagnetic
    + ✔ Plane wave
    + ✔ Field of electric/magnetic ring current
    + ✔ Field of electric/magnetic dipole
    + ✔ TE/TM spherical vector waves
    + ✔ Uniform static electric field
    + ⌛ Static charge(s)

- Acoustic
    + ✔ Plane wave
    + ✔ Monopole

##### Available scattering objects:

- Electromagnetic

    - ✔ PEC sphere
    - ⌛ PMC sphere
    - ⌛ Dielectric sphere 
    - ⌛ Multilayer dielectric sphere 
    - ⌛ Multilayer dielectric sphere with PEC core 
    - ✔ Dielectric sphere with thin impedance layer

- Acoustic

    - ✔ Sound-hard/soft sphere
    - ✔ Sound-hard/soft prolate spheroid
    - ✔ Sound-hard/soft oblate spheroid
    - ✔ Sound-hard/soft disc

##### Available quantities (where applicable):
- ✔ Far-fields
- ✔ Near-fields (electric & magnetic)
- ✔ Radar cross section (RCS)
- ⌛ Surface currents
- ✔ Scalar potentials 
- ✔ Displacement fields 

- ✔ Pressure
- ✔ Pressure traces



        
