
# the quantities obtained from the series for the scattered pressure itself; the Neumann trace
# requires the derivatives of that series and is handled separately
const AcousticQuantity = Union{Pressure,FarField,PressureTrace}


"""
    scatteredfield(sphere::Sphere, excitation::AcousticPlaneWave, quantity::AcousticQuantity; parameter::Parameter=Parameter())

Compute the pressure, the far field, or the Dirichlet trace scattered by a sound-hard or sound-soft sphere, for
an incident acoustic plane wave.
"""
function scatteredfield(sphere::Sphere, excitation::AcousticPlaneWave, quantity::AcousticQuantity; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    # no rotation of the coordinate system is required: the scattered field is a scalar
    # and rotationally symmetric about the direction of incidence

    p = progress(length(quantity.locations))

    # --- compute field
    @tasks for ind in eachindex(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, quantity.locations[ind], quantity; parameter=parameter)
        next!(p)
    end
    finish!(p)

    return F
end



"""
    scatteredfield(sphere::Sphere, excitation::AcousticPlaneWave, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace of the pressure scattered by a sound-hard or sound-soft sphere, for an incident
acoustic plane wave, employing the normals of `quantity`.
"""
function scatteredfield(
    sphere::Sphere, excitation::AcousticPlaneWave, quantity::PressureNormalGradient; parameter::Parameter=Parameter()
)

    T = typeof(excitation.frequency)
    F = zeros(Complex{T}, size(quantity.locations))

    p = progress(length(quantity.locations))

    # --- compute trace
    @tasks for ind in eachindex(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, quantity.locations[ind], quantity.normals[ind], quantity; parameter=parameter)
        next!(p)
    end
    finish!(p)

    return F
end



"""
    scatteredfield(sphere::Union{HardSphere,SoftSphere}, excitation::AcousticPlaneWave, point, quantity::AcousticQuantity; parameter::Parameter=Parameter())

Compute the pressure, the far field, or a surface trace scattered by a sound-hard or sound-soft sphere, for an
incident acoustic plane wave.

With the incident pressure ``p_\\mathrm{i} = A \\mathrm{e}^{-\\mathrm{j} k \\hat{d} ⋅ \\mathbf{r}}`` expanded as
``p_\\mathrm{i} = A \\sum_n (2n+1) (-\\mathrm{j})^n j_n(kr) P_n(\\cos\\vartheta)``, the scattered pressure reads
``p_\\mathrm{s} = A \\sum_n (2n+1) (-\\mathrm{j})^n b_n h_n^{(2)}(kr) P_n(\\cos\\vartheta)``, where ``\\vartheta`` is
measured from the direction of incidence ``\\hat{d}``. Note that, in contrast to the electromagnetic case, the
monopole term ``n = 0`` contributes.

The traces are evaluated on the surface of the sphere: only the direction of the point is taken into account,
the radial coordinate is replaced by the radius of the sphere. Hence, the locations may also be given by points
of a surface mesh, which do not lie exactly on the sphere.

The point is in Cartesian coordinates.
"""
function scatteredfield(
    sphere::Union{HardSphere,SoftSphere},
    excitation::AcousticPlaneWave,
    point,
    quantity::AcousticQuantity;
    parameter::Parameter=Parameter(),
)

    k = wavenumber(excitation)
    T = typeof(k)

    eps = parameter.relativeAccuracy

    r = norm(point)

    # the far field is determined by the direction of observation alone, whereas the pressure vanishes
    # inside the sphere
    quantity isa Pressure && r < sphere.radius && return Complex{T}(0.0)

    # the traces are evaluated on the surface of the sphere
    kr = quantity isa Trace ? k * sphere.radius : k * r

    # --- cosine of the angle enclosed by the direction of incidence and the observation direction
    cosϑ = iszero(r) ? T(1.0) : clamp(dot(excitation.direction, point) / r, -T(1.0), T(1.0))

    u = Complex{T}(0.0) # initialize
    δu = T(Inf)

    # first two values of the Legendre polynomial: P₋₁ is only required to start the recurrence
    Pn₋₁ = T(0.0)
    Pn = T(1.0)

    n = -1 # the monopole term n = 0 contributes

    try
        while δu > eps || n < 10
            n += 1

            bₙ = scatterCoeff(sphere, excitation, n)
            R = expansion(excitation, quantity, kr, n)

            Δu = (2 * n + 1) * bₙ * R * Pn

            u += Δu

            # a vanishing Legendre polynomial (e.g. Pₙ(0) for odd n) must not be mistaken for
            # convergence, hence the relative change is only updated for non-vanishing terms
            if !iszero(Δu)
                δu = iszero(u) ? T(Inf) : abs(Δu) / abs(u) # relative change
            end

            # recurrence relationship for the next Legendre polynomial
            Pn₋₁, Pn = Pn, ((2 * n + 1) * cosϑ * Pn - n * Pn₋₁) / (n + 1)
        end
    catch
        print("did not converge: n=$n\n") # if Hankel function throws overflow error
    end

    return amplitude(excitation, quantity) * u
end



"""
    scatteredfield(sphere::Union{HardSphere,SoftSphere}, excitation::AcousticPlaneWave, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace ``γ_1 p_\\mathrm{s} = \\hat{n} ⋅ ∇ p_\\mathrm{s}`` of the pressure scattered by a
sound-hard or sound-soft sphere for the given `normal`.

Decomposing the gradient in the spherical basis about the direction of incidence, where the scattered pressure
does not depend on ``φ``, yields

```math
\\hat{n} ⋅ ∇ p_\\mathrm{s} = (\\hat{n} ⋅ \\hat{r}) \\frac{∂ p_\\mathrm{s}}{∂r}
                           + (\\hat{n} ⋅ \\hat{\\vartheta}) \\frac{1}{r} \\frac{∂ p_\\mathrm{s}}{∂\\vartheta}
```

with ``\\hat{\\vartheta} = (\\cos\\vartheta \\, \\hat{r} - \\hat{d}) / \\sin\\vartheta``. Since
``∂P_n(\\cos\\vartheta) / ∂\\vartheta = -\\sin\\vartheta \\, P_n'(\\cos\\vartheta)``, the factor
``\\sin\\vartheta`` cancels, so that no special treatment of the poles is required. Differentiating the series
for the scattered pressure term by term then gives

```math
γ_1 p_\\mathrm{s} = A \\sum_n (2n+1) (-\\mathrm{j})^n b_n
    \\left[ c_r \\, k \\, h_n^{(2)\\prime}(ka) P_n(\\cos\\vartheta) - c_\\vartheta \\, h_n^{(2)}(ka) P_n'(\\cos\\vartheta) \\right]
```

where ``c_r = \\hat{n} ⋅ \\hat{r}`` and ``c_\\vartheta = (\\cos\\vartheta \\, c_r - \\hat{n} ⋅ \\hat{d}) / a``. For
the outward normal ``\\hat{n} = \\hat{r}`` the coefficient ``c_\\vartheta`` vanishes and the radial derivative
remains.

The trace is evaluated on the surface of the sphere: only the direction of the point is taken into account. The
point and the normal are in Cartesian coordinates, the latter being a unit vector.
"""
function scatteredfield(
    sphere::Union{HardSphere,SoftSphere},
    excitation::AcousticPlaneWave,
    point,
    normal,
    quantity::PressureNormalGradient;
    parameter::Parameter=Parameter(),
)

    k = wavenumber(excitation)
    T = typeof(k)

    eps = parameter.relativeAccuracy

    a  = sphere.radius # the trace is evaluated on the surface of the sphere
    ka = k * a
    s  = sqrt(π / 2 / ka)

    d = excitation.direction
    r̂ = normalize(point)

    cosϑ = clamp(dot(d, r̂), -T(1.0), T(1.0))

    # --- components of the normal along r̂ and ϑ̂, the latter including the factor 1/r
    cᵣ = dot(normal, r̂)
    cϑ = (cosϑ * cᵣ - dot(normal, d)) / a

    u  = Complex{T}(0.0) # initialize
    δu = T(Inf)

    # first two values of the Legendre polynomial and of its derivative: the entries for n = -1 are
    # only required to start the recurrences
    Pn₋₁, Pn   = T(0.0), T(1.0)
    dPn₋₁, dPn = T(0.0), T(0.0)

    n = -1 # the monopole term n = 0 contributes

    try
        while δu > eps || n < 10
            n += 1

            bₙ = scatterCoeff(sphere, excitation, n)

            h  = s * hankelh2(n + T(0.5), ka)  # spherical Hankel function
            h2 = s * hankelh2(n - T(0.5), ka)  # for derivative needed

            # Use recurrence relationship fₙ′(x) = fₙ₋₁(x) - (n + 1) / x * fₙ(x)
            dh = h2 - (n + 1) / ka * h  # derivative spherical Hankel function

            Δu = (2 * n + 1) * (-im)^n * bₙ * (cᵣ * k * dh * Pn - cϑ * h * dPn)

            u += Δu

            # a vanishing Legendre polynomial (e.g. Pₙ(0) for odd n) must not be mistaken for
            # convergence, hence the relative change is only updated for non-vanishing terms
            if !iszero(Δu)
                δu = iszero(u) ? T(Inf) : abs(Δu) / abs(u) # relative change
            end

            # recurrence relationships for the next Legendre polynomial and its derivative
            Pn₋₁, Pn, dPn₋₁, dPn = Pn,
            ((2 * n + 1) * cosϑ * Pn - n * Pn₋₁) / (n + 1), dPn,
            ((2 * n + 1) * (Pn + cosϑ * dPn) - n * dPn₋₁) / (n + 1)
        end
    catch
        print("did not converge: n=$n\n") # if Hankel function throws overflow error
    end

    return amplitude(excitation, quantity) * u
end



"""
    scatteredfield(sphere::Sphere, excitation::AcousticPlaneWave, point, quantity::AcousticQuantity; parameter::Parameter=Parameter())

Descriptive error for the field scattered by spheres for which no acoustic solution is implemented.
"""
function scatteredfield(
    sphere::Sphere, excitation::AcousticPlaneWave, point, quantity::AcousticQuantity; parameter::Parameter=Parameter()
)

    return error("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")
end



"""
    scatteredfield(sphere::Sphere, excitation::AcousticPlaneWave, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Descriptive error for the trace scattered by spheres for which no acoustic solution is implemented.
"""
function scatteredfield(
    sphere::Sphere, excitation::AcousticPlaneWave, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter()
)

    return error("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")
end



"""
    expansion(excitation::AcousticPlaneWave, quantity::Union{Pressure,PressureTrace}, kr, n::Int)

Compute the radial dependence ``(-\\mathrm{j})^n h_n^{(2)}(kr)`` of the n-th term of the series for the scattered pressure.

The Dirichlet trace shares this dependence: it is the scattered pressure evaluated at ``r = a``.
"""
function expansion(excitation::AcousticPlaneWave, quantity::Union{Pressure,PressureTrace}, kr, n::Int)

    T = typeof(excitation.frequency)

    s = sqrt(π / 2 / kr)

    return (-im)^n * s * hankelh2(n + T(0.5), kr) # spherical Hankel function
end



"""
    expansion(excitation::AcousticPlaneWave, quantity::FarField, kr, n::Int)

Compute the radial dependence of the n-th term of the series for the scattered far field.

Since ``h_n^{(2)}(kr) → \\mathrm{j}^{n+1} \\mathrm{e}^{-\\mathrm{j}kr} / (kr)`` for ``kr → ∞``, and since
``(-\\mathrm{j})^n \\mathrm{j}^{n+1} = \\mathrm{j}`` for every ``n``, the radial dependence is the same for all
terms of the series. The common factor is collected in [`amplitude`](@ref) and in the omitted
``\\mathrm{e}^{-\\mathrm{j}kr} / r``, so that unity remains.
"""
function expansion(excitation::AcousticPlaneWave, quantity::FarField, kr, n::Int)

    return Complex{typeof(excitation.frequency)}(1.0)
end



"""
    amplitude(excitation::AcousticPlaneWave, quantity::Union{Pressure,Trace})

Returns the amplitude ``A`` of the incident acoustic plane wave.
"""
amplitude(excitation::AcousticPlaneWave, quantity::Union{Pressure,Trace}) = excitation.amplitude



"""
    amplitude(excitation::AcousticPlaneWave, quantity::FarField)

Returns ``\\mathrm{j}A/k``, the factor left over by the far-field limit of the spherical Hankel functions.

As for the electromagnetic excitations, the factor ``\\mathrm{e}^{-\\mathrm{j}kr} / r`` is not included in the
far field, that is, the returned quantity is ``\\lim_{r → ∞} r \\, \\mathrm{e}^{\\mathrm{j}kr} p_\\mathrm{s}``.
"""
amplitude(excitation::AcousticPlaneWave, quantity::FarField) = im / wavenumber(excitation) * excitation.amplitude



"""
    scatterCoeff(sphere::HardSphere, excitation::AcousticPlaneWave, n::Int)

Compute the expansion coefficient ``b_n`` of the field scattered by a sound-hard sphere.

The normal velocity, and hence the radial derivative of the total pressure, vanishes on the surface, so that
``j_n'(ka) + b_n h_n^{(2)\\prime}(ka) = 0``.
"""
function scatterCoeff(sphere::HardSphere, excitation::AcousticPlaneWave, n::Int)

    T = typeof(excitation.frequency)

    ka = wavenumber(excitation) * sphere.radius

    s = sqrt(π / 2 / ka)

    j = s * besselj(n + T(0.5), ka)   # spherical Bessel function
    h = s * hankelh2(n + T(0.5), ka)  # spherical Hankel function
    j2 = s * besselj(n - T(0.5), ka)   # for derivative needed
    h2 = s * hankelh2(n - T(0.5), ka)  # for derivative needed

    # Use recurrence relationship fₙ′(x) = fₙ₋₁(x) - (n + 1) / x * fₙ(x)
    dj = j2 - (n + 1) / ka * j  # derivative spherical Bessel function
    dh = h2 - (n + 1) / ka * h  # derivative spherical Hankel function

    return -dj / dh
end



"""
    scatterCoeff(sphere::SoftSphere, excitation::AcousticPlaneWave, n::Int)

Compute the expansion coefficient ``b_n`` of the field scattered by a sound-soft sphere.

The total pressure vanishes on the surface, so that ``j_n(ka) + b_n h_n^{(2)}(ka) = 0``.
"""
function scatterCoeff(sphere::SoftSphere, excitation::AcousticPlaneWave, n::Int)

    T = typeof(excitation.frequency)

    ka = wavenumber(excitation) * sphere.radius

    s = sqrt(π / 2 / ka)

    j = s * besselj(n + T(0.5), ka)   # spherical Bessel function
    h = s * hankelh2(n + T(0.5), ka)  # spherical Hankel function

    return -j / h
end
