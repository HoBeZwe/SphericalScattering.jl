
"""
    field(excitation::AcousticMonopole, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure of an acoustic monopole, or its Dirichlet trace.

The Dirichlet trace is the pressure itself. Its locations are taken as they are given and are assumed to lie on
the surface of interest, whereas the traces of the scattered and of the total field are evaluated on the surface
of the sphere.
"""
function field(excitation::AcousticMonopole, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter(), zeroRadius=0.0)

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
    field(excitation::AcousticMonopole, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace of the pressure of an acoustic monopole for the normals of `quantity`.

The locations are taken as they are given and are assumed to lie on the surface of interest, whereas the traces
of the scattered and of the total field are evaluated on the surface of the sphere.
"""
function field(excitation::AcousticMonopole, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    F = zeros(Complex{T}, size(quantity.locations))

    # --- compute trace
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity.normals[ind], quantity; parameter=parameter)
    end

    return F
end



"""
    field(excitation::AcousticMonopole, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

Compute the pressure ``p_\\mathrm{i} = A \\mathrm{e}^{-\\mathrm{j}kR} / (4πR)`` of an acoustic monopole, which is
at the same time its Dirichlet trace ``γ_0 p_\\mathrm{i}``.

Here, ``R`` denotes the distance to the monopole, so that the pressure is the free-space Green's function of the
Helmholtz equation scaled by the amplitude ``A``. The field is singular at the position of the monopole.

The point is in Cartesian coordinates.
"""
function field(excitation::AcousticMonopole, point, quantity::Union{Pressure,PressureTrace}; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    R = norm(point - excitation.position)

    return a / (4 * π) * cis(-k * R) / R
end



"""
    field(excitation::AcousticMonopole, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

Compute the Neumann trace ``γ_1 p_\\mathrm{i} = \\hat{n} ⋅ ∇ p_\\mathrm{i}`` of the pressure of an acoustic
monopole for the given `normal`.

Since the pressure depends on the position solely via the distance ``R`` to the monopole,

```math
∇ p_\\mathrm{i} = \\frac{∂ p_\\mathrm{i}}{∂R} \\hat{R}
                = -\\frac{A}{4π} (1 + \\mathrm{j}kR) \\frac{\\mathrm{e}^{-\\mathrm{j}kR}}{R^2} \\hat{R}
```

holds, where ``\\hat{R}`` points from the monopole to the point. For ``kR ≫ 1`` this approaches
``-\\mathrm{j} k (\\hat{R} ⋅ \\hat{n}) p_\\mathrm{i}``, the local plane-wave result.

The point and the normal are in Cartesian coordinates, the latter being a unit vector.
"""
function field(excitation::AcousticMonopole, point, normal, quantity::PressureNormalGradient; parameter::Parameter=Parameter())

    a = excitation.amplitude
    k = wavenumber(excitation)

    d = point - excitation.position # vector pointing from the monopole to the point
    R = norm(d)

    return -a / (4 * π) * (1 + im * k * R) * cis(-k * R) / R^2 * dot(normal, d / R)
end
