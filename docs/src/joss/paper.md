---
title: 'SphericalScattering: A Julia Package for Electromagnetic Scattering from Spherical Objects'
tags:
  - Julia
  - dipole
  - electromagnetics
  - field theory
  - Mie scattering
  - plane wave
  - ring-current
  - sphere
  - spherical waves
  - time-harmonic
authors:
  - name: Bernd Hofmann
    orcid: 0000-0003-1435-6203
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Paula Respondek
    orcid: 0009-0005-4892-2711
    affiliation: 2
  - name: Simon B. Adrian
    orcid: 0000-0001-8008-6235
    affiliation: 2
affiliations:
 - name: Department of Electrical Engineering, School of Computation, Information and Technology, Technical University of Munich, 80290 Munich, Germany
   index: 1
 - name: Fakultät für Informatik und Elektrotechnik, Universität Rostock, 18059 Rostock, Germany
   index: 2
software_repository_url:
 - https://github.com/HoBeZwe/SphericalScattering.jl
date: 13 June 2023
bibliography: paper.bib
---

# Summary   

When electromagnetic fields are impinging on objects of various kinds, determining the scattered field as a solution to Maxwell's equations is crucial for many applications.
For example, when monitoring the position of an airplane by a radar, the scattering behavior of the airplane plays a pivotal role and, thus, needs to be studied.
Analytical approaches, however, to characterize such scattering behavior are rarely known.
Some of the few exceptions where at least semi-analytical descriptions are available are metallic or dielectric spherical objects excited by time-harmonic or static fields [@ruckRadarCrossSection1970;@jinTheoryComputationElectromagnetic2015].
While there are some applications where canonical scattering problems are the study subject of interest, in other areas solutions to the scattering from spherical objects rather serve as a means to verify the correctness of more involved numerical techniques.
The latter allow to analyze the scattering from real-world objects, for instance, via finite element or integral equation methods [@raoElectromagneticScatteringSurfaces1982;@harringtonFieldComputationMoment1993;@jinTheoryComputationElectromagnetic2015;@adrianElectromagneticIntegralEquations2021].


# Statement of need

`SphericalScattering` is a Julia package [@bezanson2017julia] providing semi-analytical solutions to the scattering of time-harmonic as well as static electromagnetic fields from spherical objects (including the Mie solutions for plane wave excitations).
To this end, series expansions are evaluated with special care to obtain accurate solutions down to the static limit.
The series expansions are based on expressing the incident and scattered fields in terms of spherical wave functions such that the boundary conditions can be enforced at interfaces of different materials yielding the expansion coefficients of the spherical wave functions of the scattered field [@jinTheoryComputationElectromagnetic2015;@ruckRadarCrossSection1970]. 

Other available implementations have a different focus, that is, specific 2D scenarios are addressed [@Blankrot2018], T-matrices are employed for general shaped objects [@art_gower_2018_1213225;@schebarchov2021;@miepy2022;@Egel2017], ensemble averaged waves are obtained [@fidgit2020], spontaneous decay rates of a dipole are studied [@Rasskazov20], light scattering is considered employing only plane waves as excitations [@Prahl_miepython_Pure_python;@konstantin_ladutenko_2017_248729;@matscat2023;@matmie2023;@piemiecoated2016;@miescatwu2023;@cppmie2022], or only far-field quantities are computed. 

In contrast, in `SphericalScattering` a variety of excitations is available, that is,

- plane waves,
- fields of electric/magnetic ring currents,
- fields of electric/magnetic dipoles,
- transverse electric (TE) and transverse magnetic (TM) spherical vector waves, and
- uniform static electric fields,

where several parameters including the orientation, direction, or polarization of the sources can be set by the user and are not predefined.
The scattered far- and near-fields are then obtained following [@jinTheoryComputationElectromagnetic2015;@ruckRadarCrossSection1970;@sihvolaTransmissionLineAnalogy1988;@hansen1988;@jackson1999;@jones1995] for

- perfectly electrically conducting (PEC) spheres and
- dielectric spheres

all via a unified interface.
In consequence, `SphericalScattering` facilitates a reproducible and comparable verification of approaches to solve electromagnetic scattering problems.
For this purpose, it has already been employed in scientific puplications [@hofmannLowFrequencyStableDiscretization2021;@hofmannLowFrequencyStabilizedElectricField2022;@hofmannEfficientCombinationScalarPotential2022a;@hofmannExcitationAwareSelfAdaptiveFrequency2023;@hofmannLowFrequencyStableExcitation2023;@hofmannInvestigationsLowFrequencyStability2023;@hofmannSelfAdaptiveFrequencyNormalization2023].


# Acknowledgments
Paula Respondek and Simon B. Adrian were funded by the Deutsche Forschungsgemeinschaft (DFG; German Research Foundation) under Grant SFB 1270/2-299150580.



# References
