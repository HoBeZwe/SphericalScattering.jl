
"""
    symmetryAxis(excitation::AcousticMonopole)

Returns the direction towards the monopole, about which the scattered field is rotationally symmetric.
"""
symmetryAxis(excitation::AcousticMonopole) = normalize(excitation.position)



"""
    incidentCoeff(excitation::AcousticMonopole, n::Int)

Compute the coefficient ``e_n = -\\mathrm{j} k \\, h_n^{(2)}(k R_0) / (4π)`` of the n-th term of the incident
expansion, where ``R_0`` denotes the distance of the monopole from the center of the sphere.

It follows from the addition theorem of the free-space Green's function,

```math
\\frac{\\mathrm{e}^{-\\mathrm{j}k|\\mathbf{r} - \\mathbf{r}_0|}}{4π|\\mathbf{r} - \\mathbf{r}_0|}
    = \\frac{-\\mathrm{j}k}{4π} \\sum_n (2n+1) j_n(k r_<) h_n^{(2)}(k r_>) P_n(\\cos\\vartheta)
```

with ``r_< = \\min(r, R_0)``, ``r_> = \\max(r, R_0)`` and ``\\vartheta`` measured from the monopole. Since the
monopole lies outside the sphere, ``r_> = R_0`` holds at its surface, which is where the boundary condition
determines the scattering coefficients.
"""
function incidentCoeff(excitation::AcousticMonopole, n::Int)

    T = typeof(excitation.frequency)

    k = wavenumber(excitation)
    kR₀ = k * norm(excitation.position)

    s = sqrt(π / 2 / kR₀)

    return -im * k / (4 * π) * s * hankelh2(n + T(0.5), kR₀) # spherical Hankel function
end



"""
    checkExcitation(sphere::Sphere, excitation::AcousticMonopole)

Ensure that the monopole lies outside the sphere, as the expansion of the incident field assumes.
"""
function checkExcitation(sphere::Sphere, excitation::AcousticMonopole)

    norm(excitation.position) > sphere.radius ||
        error("The monopole has to be located outside the sphere: its distance from the center is smaller than the radius.")

    return nothing
end
