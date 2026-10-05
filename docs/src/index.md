
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

To install julia, you can follow [these instructions](https://docs.julialang.org/en/v1/manual/getting-started/).
Installing SphericalScattering is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add SphericalScattering 
```


---
## Feature Overview

The following aspects are implemented (✔) and planned (⌛):

```@raw html
<style>
  #documenter .content table.feature-table {
    display: table; /* Documenter makes tables blocks, whose borders would span the whole page */
    width: auto;
    border-collapse: collapse;
    border-top: 2px solid currentColor;
    border-bottom: 2px solid currentColor;
  }
  #documenter .content table.feature-table th,
  #documenter .content table.feature-table td { vertical-align: middle; padding: 0.35em 0.6em; }
  #documenter .content table.feature-table thead th { text-align: center; }
  #documenter .content table.feature-table thead th.ft-left { text-align: left; }
  #documenter .content table.feature-table tr.ft-group th {
    text-align: left;
    font-weight: 600;
    letter-spacing: 0.02em;
    background: rgba(127, 127, 127, 0.12);
    border-top: 2px solid rgba(127, 127, 127, 0.35);
  }
  #documenter .content table.feature-table td.ft-object { padding-left: 1.1em; }
  #documenter .content table.feature-table td.ft-check { text-align: center; }
  #documenter .content table.feature-table .ft-note { opacity: 0.7; font-size: 0.85em; }
</style>
<div style="overflow-x: auto;">
<table class="feature-table">
  <thead>
    <tr>
      <th class="ft-left" rowspan="2">object</th>
      <th class="ft-left" rowspan="2">boundary</th>
      <th colspan="8">electromagnetic excitation</th>
    </tr>
    <tr>
      <th>plane wave</th>
      <th>el. ring current</th>
      <th>mag. ring current</th>
      <th>el. dipole</th>
      <th>mag. dipole</th>
      <th>TE/TM modes</th>
      <th>uniform static field</th>
      <th>static charge(s)</th>
    </tr>
  </thead>
  <tbody>
    <tr class="ft-group"><th colspan="10">Spheres</th></tr>
    <tr><td class="ft-object" rowspan="6"> </td><td>PEC</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">✔</td><td class="ft-check">⌛</td></tr>
    <tr><td>PMC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td>dielectric</td><td class="ft-check">✔</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">✔</td><td class="ft-check">⌛</td></tr>
    <tr><td>multilayer dielectric</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">✔</td><td class="ft-check">⌛</td></tr>
    <tr><td>multilayer dielectric with PEC core</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">✔</td><td class="ft-check">⌛</td></tr>
    <tr><td>dielectric with thin impedance layer</td><td class="ft-check">➖</td><td class="ft-check">➖</td><td class="ft-check">➖</td><td class="ft-check">➖</td><td class="ft-check">➖</td><td class="ft-check">➖</td><td class="ft-check">✔</td><td class="ft-check">➖</td></tr>

    <tr class="ft-group"><th colspan="10">Spheroids</th></tr>
    <tr><td class="ft-object" rowspan="2">prolate</td><td>PEC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td>PMC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td class="ft-object" rowspan="2">oblate</td><td>PEC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td>PMC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td class="ft-object" rowspan="2">disc<br/><span class="ft-note">flat oblate spheroid</span></td><td>PEC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
    <tr><td>PMC</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td><td class="ft-check">⌛</td></tr>
  </tbody>
</table>
</div>
```
---

```@raw html
<div style="overflow-x: auto;">
<table class="feature-table">
  <thead>
    <tr>
      <th class="ft-left" rowspan="2">object</th>
      <th class="ft-left" rowspan="2">boundary</th>
      <th colspan="2">acoustic excitation</th>
    </tr>
    <tr>
      <th>plane wave</th>
      <th>monopole</th>
    </tr>
  </thead>
  <tbody>
    <tr class="ft-group"><th colspan="4">Spheres</th></tr>
    <tr><td class="ft-object" rowspan="2"> </td><td>sound-hard</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td>sound-soft</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>

    <tr class="ft-group"><th colspan="4">Spheroids</th></tr>
    <tr><td class="ft-object" rowspan="2">prolate </td><td>sound-hard</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td>sound-soft</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td class="ft-object" rowspan="2">oblate </td><td>sound-hard</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td>sound-soft</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td class="ft-object" rowspan="2">disc<br/><span class="ft-note">flat oblate spheroid</span></td><td>sound-hard</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
    <tr><td>sound-soft</td><td class="ft-check">✔</td><td class="ft-check">✔</td></tr>
  </tbody>
</table>
</div>
```


---
##### Available quantities

The quantities which can be computed for these setups, among them near fields, far fields, potentials, and surface traces, are listed under [Quantities](@ref quantitiesConcept). The surface currents of the electromagnetic scatterers are planned (⌛).
