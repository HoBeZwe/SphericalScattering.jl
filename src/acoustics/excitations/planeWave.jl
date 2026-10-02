


struct AcousticPlaneWave{T,R,C} <: AcousticExcitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    direction::SVector{3,R}

    # inner constructor: normalize direction
    function AcousticPlaneWave(embedding::Medium{C}, frequency::R, amplitude::T, direction::SVector{3,R}) where {T,R,C}

        dir_normalized = normalize(direction)

        new{T,R,C}(embedding, frequency, amplitude, dir_normalized)
    end
end



"""
    field(excitation::AcousticPlaneWave, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure of an acoustic plane wave, or its Dirichlet trace.

The Dirichlet trace is the pressure itself. Its locations are taken as they are given and are assumed to lie on
the surface of interest, whereas the traces of the scattered and of the total field are evaluated on the surface
of the sphere.
"""
function field(
    excitation::AcousticPlaneWave, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter(), zeroRadius=0.0
)

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        norm(point) < zeroRadius && continue
        F[ind] = field(excitation, point, quantity; parameter=parameter)
    end

    return F
end



"""
    field(excitation::AcousticPlaneWave, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace of the pressure of an acoustic plane wave for the normals of `quantity`.

The locations are taken as they are given and are assumed to lie on the surface of interest, whereas the traces
of the scattered and of the total field are evaluated on the surface of the sphere.
"""
function field(excitation::AcousticPlaneWave, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute trace
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity.normals[ind], quantity; parameter=parameter)
    end

    return F
end



"""
    field(excitation::AcousticPlaneWave, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure ``p_\\mathrm{i} = A \\mathrm{e}^{-\\mathrm{j} k \\hat{d} ⋅ \\mathbf{r}}`` of an acoustic plane
wave, which is at the same time its Dirichlet trace ``γ_0 p_\\mathrm{i}``.

The point is in Cartesian coordinates.
"""
function field(excitation::AcousticPlaneWave, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    d = excitation.direction

    return a * cis(-k * dot(d, point))
end



"""
    field(excitation::AcousticPlaneWave, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace ``γ_1 p_\\mathrm{i} = \\hat{n} ⋅ ∇ p_\\mathrm{i}`` of the pressure of an acoustic plane
wave for the given `normal`.

Since ``∇ p_\\mathrm{i} = -\\mathrm{j} k \\hat{d} \\, p_\\mathrm{i}``, the trace evaluates to
``-\\mathrm{j} k (\\hat{d} ⋅ \\hat{n}) \\, p_\\mathrm{i}``.

The point and the normal are in Cartesian coordinates, the latter being a unit vector.
"""
function field(excitation::AcousticPlaneWave, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    d = excitation.direction

    return -im * k * dot(d, normal) * a * cis(-k * dot(d, point))
end



"""
    field(excitation::AcousticPlaneWave, quantity::FarField; parameter::Parameter=Parameter())

Throw error since the far-field of a plane wave is not defined.
"""
function field(excitation::AcousticPlaneWave, quantity::FarField; parameter::Parameter=Parameter(), zeroRadius=0.0)

    return error("The far-field of a plane wave is not defined.")
end



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
