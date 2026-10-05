
# The radial functions of SpheroidalWaves.jl come in four kinds: kind 1 and 2 are the two real solutions
# R⁽¹⁾ and R⁽²⁾, whereas kind 3 and 4 are their complex combinations R⁽¹⁾ ± j R⁽²⁾. This package employs
# the time convention e^{jωt}, for which the outgoing solution is the analogue of hₙ⁽²⁾ = jₙ - j yₙ, that
# is, kind 4. Note that the spheroidal literature calls the outgoing function R⁽³⁾, which corresponds to
# the opposite time convention: picking kind 3 would yield an incoming wave. The convention concerns the
# time dependence alone and holds for both shapes, as their radiation conditions confirm.
const regularKind = 1
const outgoingKind = 4


"""
    spheroidalRadial(shape::Symbol, m::Int, n::Int, c, ξ; outgoing::Bool=false)

Compute the spheroidal radial function of the given `shape`, `:oblate` or `:prolate`, of order `m` and degree
`n` and its derivative with respect to ``ξ`` at the spheroidal parameter `c` and the radial coordinate `ξ`.

Returned is a named tuple `(value, derivative)`. For `outgoing=false` the regular function ``R^{(1)}_{mn}`` is
evaluated, which is the counterpart of ``j_n``; for `outgoing=true` the outgoing one is evaluated, which is the
counterpart of ``h_n^{(2)}``.
"""
function spheroidalRadial(shape::Symbol, m::Int, n::Int, c, ξ; outgoing::Bool=false)

    kind = outgoing ? outgoingKind : regularKind

    r = rmn(m, n, c, [ξ]; spheroid=shape, kind=kind)

    return (value=only(vec(r.value)), derivative=only(vec(r.derivative)))
end


"""
    spheroidalAngular(shape::Symbol, m::Int, n::Int, c, η)

Compute the spheroidal angular function of the given `shape`, `:oblate` or `:prolate`, of order `m` and degree
`n` and its derivative with respect to ``η`` at the spheroidal parameter `c` and the angular coordinate `η`.

Returned is a named tuple `(value, derivative)`. The Meixner-Schäfke normalization is employed, for which the
angular function reduces to the associated Legendre function ``P_n^m`` as ``c → 0``, so that the limit of a
sphere connects directly to the series of the spherical scatterers.
"""
function spheroidalAngular(shape::Symbol, m::Int, n::Int, c, η)

    s = smn(m, n, c, [η]; spheroid=shape, normalize=false)

    return (value=only(vec(s.value)), derivative=only(vec(s.derivative)))
end


"""
    spheroidalParameter(sphere::Spheroid, excitation::Excitation)

Returns the spheroidal parameter ``c = k f``, the product of the wavenumber and the semifocal distance, which
takes the role of ``ka`` for a sphere.
"""
spheroidalParameter(sphere::Spheroid, excitation::Excitation) = wavenumber(excitation) * sphere.semifocal


"""
    angularNorm(m::Int, n::Int)

Returns the norm ``N_{mn} = \\int_{-1}^{1} S_{mn}^2(c, η) \\, \\mathrm{d}η = \\cfrac{2}{2n + 1} \\cfrac{(n + m)!}{(n - m)!}``
of the angular functions in the Meixner-Schäfke normalization.

It does not depend on ``c``: the normalization is defined such that the norm is that of the associated Legendre
functions, to which the angular functions reduce as ``c → 0``. The ratio of the factorials is formed as a product,
which is exact up to ``2m`` roundings.
"""
function angularNorm(m::Int, n::Int)

    ratio = 1.0
    for k in (n - m + 1):(n + m)
        ratio *= k
    end

    return 2 / (2n + 1) * ratio
end
