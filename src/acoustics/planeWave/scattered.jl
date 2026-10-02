
"""
    symmetryAxis(excitation::AcousticPlaneWave)

Returns the direction of incidence, about which the scattered field is rotationally symmetric.
"""
symmetryAxis(excitation::AcousticPlaneWave) = excitation.direction



"""
    incidentCoeff(excitation::AcousticPlaneWave, n::Int)

Compute the coefficient ``e_n = (-\\mathrm{j})^n`` of the n-th term of the incident expansion.

It follows from ``\\mathrm{e}^{-\\mathrm{j} k \\hat{d} ⋅ \\mathbf{r}}
= \\sum_n (2n+1) (-\\mathrm{j})^n j_n(kr) P_n(\\cos\\vartheta)``, the plane-wave expansion with ``\\vartheta``
measured from the direction of incidence ``\\hat{d}``.
"""
incidentCoeff(excitation::AcousticPlaneWave, n::Int) = (-im)^n
