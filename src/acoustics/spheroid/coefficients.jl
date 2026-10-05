
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

    R₁ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀)
    Rₒ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀; outgoing=true)

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

    R₁ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀)
    Rₒ = spheroidalRadial(shape(sphere), m, n, c, sphere.ξ₀; outgoing=true)

    return -R₁.value / Rₒ.value
end
