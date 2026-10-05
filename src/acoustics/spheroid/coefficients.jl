
"""
    scatterCoeff(sphere::Spheroid, excitation::AcousticExcitation, m::Int, n::Int)

Compute the expansion coefficient ``b_{mn}`` of the field scattered by a spheroid.

Since the surface ``ξ = ξ_0`` is a coordinate surface and the angular functions form a complete orthogonal set
on it, the boundary condition decouples the modes, see [`boundaryRatio`](@ref). As for a sphere, the coefficient
does not depend on the excitation.
"""
function scatterCoeff(sphere::Spheroid, excitation::AcousticExcitation, m::Int, n::Int)

    c = spheroidalParameter(sphere, excitation)

    R₁ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀)
    Rₒ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀; outgoing=true)

    return boundaryRatio(sphere, R₁.value, R₁.derivative, Rₒ.value, Rₒ.derivative)
end


"""
    boundaryRatio(sphere::Spheroid{SoundHard}, R₁, R₁′, Rₒ, Rₒ′)

Compute the scattering coefficient ``b_{mn}`` of a sound-hard spheroid from the regular and the outgoing radial
function and their derivatives with respect to ``ξ`` on the surface.

The normal derivative of the total pressure vanishes, so that
``R^{(1)\\prime}_{mn}(c, ξ_0) + b_{mn} R^{(\\mathrm{out})\\prime}_{mn}(c, ξ_0) = 0`` holds.

For a disc, ``ξ_0 = 0``, the coefficient vanishes identically for even ``n - m``: at the degenerate surface the
radial functions split by parity, and the Neumann problem is carried by the modes with odd ``n - m`` alone.
"""
boundaryRatio(sphere::Spheroid{SoundHard}, R₁, R₁′, Rₒ, Rₒ′) = -R₁′ / Rₒ′

"""
    boundaryRatio(sphere::Spheroid{SoundSoft}, R₁, R₁′, Rₒ, Rₒ′)

Compute the scattering coefficient ``b_{mn}`` of a sound-soft spheroid from the regular and the outgoing radial
function on the surface.

The total pressure vanishes, so that ``R^{(1)}_{mn}(c, ξ_0) + b_{mn} R^{(\\mathrm{out})}_{mn}(c, ξ_0) = 0`` holds.

For a disc, ``ξ_0 = 0``, the coefficient vanishes identically for odd ``n - m``: the Dirichlet problem is carried
by the modes with even ``n - m`` alone.
"""
boundaryRatio(sphere::Spheroid{SoundSoft}, R₁, R₁′, Rₒ, Rₒ′) = -R₁ / Rₒ



"""
    incidentCoefficients(sphere::Spheroid, excitation::AcousticExcitation, N::Int)

Compute the coefficients ``A_{mn}`` of the expansion of the incident field
``p_\\mathrm{i} = \\sum_{mn} A_{mn} R^{(1)}_{mn}(c, ξ) S_{mn}(c, η) \\mathrm{e}^{\\mathrm{j}mφ}`` in the frame of the
spheroid, for all orders ``|m| ≤ N`` and degrees ``|m| ≤ n ≤ N``, indexed as `[m + N + 1, n + 1]`.

The plane wave and the monopole have analytic expansions, see their methods. For any other excitation, whose
expansion is not known, the coefficients are obtained by projecting its field, see [`projectedCoefficients`](@ref).
"""
incidentCoefficients(sphere::Spheroid, excitation::AcousticExcitation, N::Int) = projectedCoefficients(sphere, excitation, N)


"""
    projectedCoefficients(sphere::Spheroid, excitation::AcousticExcitation, N::Int; ξ=nothing, nη=nothing, nφ=nothing, eps=1e-12)

Compute the coefficients of the incident expansion, see [`incidentCoefficients`](@ref), by projecting the incident
field onto the angular functions on a surface ``ξ = \\mathrm{const}``.

The angular functions of equal order are orthogonal there with unit weight, and the exponentials are orthogonal in
``φ``, so that

```math
A_{mn} = \\cfrac{1}{R^{(1)}_{mn}(c, ξ) \\, N_{mn}} \\int_{-1}^{1} \\left[ \\cfrac{1}{2π} \\int_0^{2π}
         p_\\mathrm{i} \\, \\mathrm{e}^{-\\mathrm{j}mφ} \\, \\mathrm{d}φ \\right] S_{mn}(c, η) \\, \\mathrm{d}η
```

holds with ``N_{mn} = \\int_{-1}^{1} S_{mn}^2 \\, \\mathrm{d}η``. The projection requires nothing but the incident
field, which makes it the fallback for excitations without an analytic expansion; for the plane wave and the
monopole it serves as an independent check of the analytic coefficients. The integrals are evaluated by
Gauss-Legendre quadrature with `nη` nodes in ``η`` and by `nφ` equispaced nodes in ``φ``; orders whose azimuthal
projection is below the relative accuracy `eps` are left zero.

The surface of the projection has to lie between the scatterer and the source of the incident field: beyond the
source, see [`sourceCoordinate`](@ref), the expansion in the regular wave functions does not hold. By default it
is placed at [`projectionCoordinate`](@ref), or halfway between the scatterer and the source if that is closer; a
given `ξ` reaching the source is rejected.

!!! note
    The coefficients of high degree are accurate only as far as they matter on the surface of the projection: the
    error of the quadrature is amplified by ``1 / R^{(1)}_{mn}(c, ξ)``, which is large for degrees beyond ``c ξ``.
    The field on and inside that surface is unaffected, since it involves the product ``A_{mn} R^{(1)}_{mn}``.
"""
function projectedCoefficients(sphere::Spheroid, excitation::AcousticExcitation, N::Int; ξ=nothing, nη=nothing, nφ=nothing, eps=1e-12)

    T = typeof(excitation.frequency)

    c = spheroidalParameter(sphere, excitation)

    # --- the surface of the projection, which has to lie between the scatterer and the source of the incident
    #     field: beyond the source the expansion in the regular wave functions does not hold
    ξs = sourceCoordinate(sphere, excitation)
    ξp = isnothing(ξ) ? T(min(projectionCoordinate(sphere), (sphere.ξ₀ + ξs) / 2)) : T(ξ)

    ξp < ξs || error(
        "The surface of the projection, ξ = $ξp, has to lie closer to the scatterer than the source of the incident field at ξ = $ξs.",
    )

    nη = isnothing(nη) ? 2 * N + 16 : nη
    nφ = isnothing(nφ) ? 4 * N + 8 : nφ

    ηNodes, ηWeights = gausslegendre(nη)
    φNodes = T(2π) .* (0:(nφ - 1)) ./ nφ

    R = frame(sphere)

    # --- the incident field on the surface of the projection, in the global frame
    quantity = Pressure(nothing) # the single-point methods take the quantity for dispatch only

    pᵢ = zeros(Complex{T}, nη, nφ)
    for (i, η) in enumerate(ηNodes), (j, φ) in enumerate(φNodes)
        pᵢ[i, j] = field(excitation, R * cartesianCoordinates(sphere, SVector(ξp, η, φ)), quantity)
    end

    # --- the azimuthal projections; orders below the relative accuracy are not projected further
    g = [(pᵢ * cis.(-m .* φNodes)) ./ nφ for m in (-N):N]
    amplitude = [maximum(abs, gₘ) for gₘ in g]
    significant = amplitude .> eps * maximum(amplitude)

    A = zeros(Complex{T}, 2 * N + 1, N + 1)

    for m in (-N):N

        significant[m + N + 1] || continue

        mAbs = abs(m)

        gₘ = g[m + N + 1]

        S = smn(mAbs, mAbs:N, c, ηNodes; spheroid=shape(sphere), normalize=false).value
        R₁ = rmn(mAbs, mAbs:N, c, [ξp]; spheroid=shape(sphere), kind=regularKind).value

        for (k, n) in enumerate(mAbs:N)
            Sₙ = view(S, :, k)

            A[m + N + 1, n + 1] = sum(ηWeights .* gₘ .* Sₙ) / (R₁[1, k] * sum(ηWeights .* Sₙ .^ 2))
        end
    end

    return A
end


"""
    sourceCoordinate(sphere::Spheroid, excitation::AcousticExcitation)

Returns the radial coordinate of the source of the incident field, beyond which its expansion in the regular
spheroidal wave functions does not hold: that of the position of a monopole, and infinity for a plane wave.
"""
sourceCoordinate(sphere::Spheroid, excitation::AcousticPlaneWave) = Inf

sourceCoordinate(sphere::Spheroid, excitation::AcousticMonopole) =
    spheroidalCoordinates(sphere, frame(sphere)' * excitation.position)[1]
