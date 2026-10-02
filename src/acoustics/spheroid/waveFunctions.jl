
# The radial functions of SpheroidalWaves.jl come in four kinds: kind 1 and 2 are the two real solutions
# R⁽¹⁾ and R⁽²⁾, whereas kind 3 and 4 are their complex combinations R⁽¹⁾ ± j R⁽²⁾. This package employs
# the time convention e^{jωt}, for which the outgoing solution is the analogue of hₙ⁽²⁾ = jₙ - j yₙ, that
# is, kind 4. Note that the spheroidal literature calls the outgoing function R⁽³⁾, which corresponds to
# the opposite time convention: picking kind 3 would yield an incoming wave.
const oblateRegularKind = 1
const oblateOutgoingKind = 4


"""
    oblateRadial(m::Int, n::Int, c, ξ; outgoing::Bool=false)

Compute the oblate spheroidal radial function of order `m` and degree `n` and its derivative with respect to
``ξ`` at the spheroidal parameter `c` and the radial coordinate `ξ`.

Returned is a named tuple `(value, derivative)`. For `outgoing=false` the regular function ``R^{(1)}_{mn}`` is
evaluated, which is the counterpart of ``j_n``; for `outgoing=true` the outgoing one is evaluated, which is the
counterpart of ``h_n^{(2)}``.
"""
function oblateRadial(m::Int, n::Int, c, ξ; outgoing::Bool=false)

    kind = outgoing ? oblateOutgoingKind : oblateRegularKind

    r = rmn(m, n, c, [ξ]; spheroid=:oblate, kind=kind)

    return (value=only(vec(r.value)), derivative=only(vec(r.derivative)))
end


"""
    oblateAngular(m::Int, n::Int, c, η)

Compute the oblate spheroidal angular function of order `m` and degree `n` and its derivative with respect to
``η`` at the spheroidal parameter `c` and the angular coordinate `η`.

Returned is a named tuple `(value, derivative)`. The Meixner-Schäfke normalization is employed, for which the
angular function reduces to the associated Legendre function ``P_n^m`` as ``c → 0``, so that the limit of a
sphere connects directly to the series of the spherical scatterers.
"""
function oblateAngular(m::Int, n::Int, c, η)

    s = smn(m, n, c, [η]; spheroid=:oblate, normalize=false)

    return (value=only(vec(s.value)), derivative=only(vec(s.derivative)))
end


"""
    spheroidalParameter(sphere::Spheroid, excitation::AcousticExcitation)

Returns the spheroidal parameter ``c = k f``, the product of the wavenumber and the semifocal distance, which
takes the role of ``ka`` for a sphere.
"""
spheroidalParameter(sphere::Spheroid, excitation::AcousticExcitation) = wavenumber(excitation) * sphere.semifocal


"""
    scatterCoeff(sphere::Spheroid{SoundHard}, excitation::AcousticExcitation, m::Int, n::Int)

Compute the expansion coefficient ``b_{mn}`` of the field scattered by a sound-hard spheroid.

Since the surface ``ξ = ξ_0`` is a coordinate surface and the angular functions form a complete orthogonal set
on it, the boundary condition decouples the modes, so that
``R^{(1)\\prime}_{mn}(c, ξ_0) + b_{mn} R^{(\\mathrm{out})\\prime}_{mn}(c, ξ_0) = 0`` holds, the prime denoting
the derivative with respect to ``ξ``. As for a sphere, the coefficient does not depend on the excitation.

For a disc, ``ξ_0 = 0``, the coefficient vanishes identically for even ``n - m``: at the degenerate surface the
radial functions split by parity, and the Neumann problem is carried by the modes with odd ``n - m`` alone.
"""
function scatterCoeff(sphere::Spheroid{SoundHard}, excitation::AcousticExcitation, m::Int, n::Int)

    c = spheroidalParameter(sphere, excitation)

    R₁ = oblateRadial(m, n, c, sphere.ξ₀)
    Rₒ = oblateRadial(m, n, c, sphere.ξ₀; outgoing=true)

    return -R₁.derivative / Rₒ.derivative
end


"""
    scatterCoeff(sphere::Spheroid{SoundSoft}, excitation::AcousticExcitation, m::Int, n::Int)

Compute the expansion coefficient ``b_{mn}`` of the field scattered by a sound-soft spheroid.

The total pressure vanishes on the surface, so that
``R^{(1)}_{mn}(c, ξ_0) + b_{mn} R^{(\\mathrm{out})}_{mn}(c, ξ_0) = 0`` holds. As for a sphere, the coefficient
does not depend on the excitation.

For a disc, ``ξ_0 = 0``, the coefficient vanishes identically for odd ``n - m``: the Dirichlet problem is
carried by the modes with even ``n - m`` alone.
"""
function scatterCoeff(sphere::Spheroid{SoundSoft}, excitation::AcousticExcitation, m::Int, n::Int)

    c = spheroidalParameter(sphere, excitation)

    R₁ = oblateRadial(m, n, c, sphere.ξ₀)
    Rₒ = oblateRadial(m, n, c, sphere.ξ₀; outgoing=true)

    return -R₁.value / Rₒ.value
end
