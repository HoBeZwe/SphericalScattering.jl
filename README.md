
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo_Scat_README_white.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo_Scat_README.svg" height="190">
  <img alt="" src="" height="190">
</picture>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hobezwe.github.io/SphericalScattering.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/HoBeZwe/SphericalScattering.jl/blob/master/LICENSE)
[![CI](https://github.com/HoBeZwe/SphericalScattering.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HoBeZwe/SphericalScattering.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/HoBeZwe/SphericalScattering.jl/branch/master/graph/badge.svg?token=4F9NUNRC1K)](https://codecov.io/gh/HoBeZwe/SphericalScattering.jl)
[![DOI](https://zenodo.org/badge/375493054.svg)](https://zenodo.org/badge/latestdoi/375493054)

## Introduction

This package provides semi-analytical solutions to the scattering of time harmonic and static electromagnetic fields from spherical objects. 
To this end, series expansions are evaluated. Special care is taken to obtain accurate solutions down to the static limit.



##### Available incident fields:
- :heavy_check_mark: Plane wave
- :heavy_check_mark: Field of electric/magnetic ring current
- :heavy_check_mark: Field of electric/magnetic dipole
- :heavy_check_mark: TE/TM spherical vector waves
- :heavy_check_mark: Uniform static electric field

##### Available scattering objects:
- :heavy_check_mark: PEC sphere
- :hourglass_flowing_sand: PMC sphere
- :hourglass_flowing_sand: Dielectric sphere (so far only for uniform static field and the plane-wave)
- :hourglass_flowing_sand: Multilayer dielectric sphere (so far only for uniform static field)
- :hourglass_flowing_sand: Multilayer dielectric sphere with PEC core (so far only for uniform static field)


## Citation

Please cite this package following the information on [Zenodo](https://zenodo.org/badge/latestdoi/375493054).


## Documentation

Here you can find the [documentation](https://hobezwe.github.io/SphericalScattering.jl/dev/).