
"""
    patternPoints(;r=1.0, resolution=5)

Convenience function returning points (in Cartesian and spherical coordinates) on a spherical grid at a distance 'r' and a resolution in degrees. 
"""
function sphericalGridPoints(; r=1.0, resolution=5, center=SVector(0.0, 0.0, 0.0))

    ϑ = range(0.0; stop=π, length=round(Int, 180 / resolution))  # default: 5° steps
    ϕ = range(0.0; stop=2π, length=round(Int, 360 / resolution)) # default: 5° steps 

    points_cart = [SVector(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) + center for θ in ϑ, φ in ϕ]
    points_sph  = [SVector(r, θ, φ) + center for θ in ϑ, φ in ϕ]

    return points_cart, points_sph
end


"""
    phiCutPoints(ϕ; r=1.0, resolution=1)

Convenience function returning points (in Cartesian coordinates) on a phi cut at a distance 'r' and a resolution in degrees. 
"""
function phiCutPoints(ϕ; r=1.0, resolution=1)

    ϑ = range(-π; stop=π, length=round(Int, 360 / resolution))  # default: 1° steps

    points_cart = [SVector(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]

    return points_cart, ϑ
end


"""
    thetaCutPoints(ϑ; r=1.0, resolution=1)

Convenience function returning points (in Cartesian coordinates) on a theta cut at a distance 'r' and a resolution in degrees. 
"""
function thetaCutPoints(ϑ; r=1.0, resolution=1)

    ϕ = range(0; stop=2π, length=round(Int, 360 / resolution))  # default: 1° steps

    points_cart = [SVector(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]

    return points_cart, ϕ
end


function surfacePoints end

# the docstrings of `surfacePoints` and `plotsurface` are attached to signatures rather than to the functions: a
# docstring of the function itself would be the very object of a `@docs` block naming the function, which the
# `@autodocs` of the API reference would then include a second time

"""
    surfacePoints(sphere::Scatterer; resolution=5)

Points on the surface of the scatterer on a grid suited for plotting surface densities, see [`plotsurface`](@ref),
with a resolution in degrees.

The points are arranged as a matrix: the rows run from the pole in the direction of the axis to the opposite one,
the columns around the axis, the last column repeating the first, so that the surface is closed. They are given in
Cartesian coordinates, so that any quantity can be computed at them by the general API, e.g.,
`field(sp, ex, PressureTrace(points))`.

For a spheroid the rows are equispaced in the polar angle ``θ`` with ``η = \\cos θ``, so that the tips of a prolate
spheroid are resolved as well as its equator. For a disc the rows with ``η > 0`` lie on the face whose normal is the
axis, those with ``η < 0`` on the opposite one; they are displaced from the disc by a tiny fraction of its radius
along the axis, which determines the face (exactly on the disc both faces share their Cartesian coordinates). The rim
itself is omitted, as the Neumann trace is singular there; the outermost rows lie half a step before it.
"""
surfacePoints(sphere::Scatterer)

function surfacePoints(sphere::Sphere; resolution=5)

    θ = range(0.0, π; length=round(Int, 180 / resolution) + 1)
    φ = range(0.0, 2π; length=round(Int, 360 / resolution) + 1)

    r = sphere.radius

    return [SVector(r * sin(t) * cos(p), r * sin(t) * sin(p), r * cos(t)) for t in θ, p in φ]
end

function surfacePoints(sphere::Spheroid; resolution=5)

    T = typeof(sphere.semifocal)

    n = round(Int, 180 / resolution)
    Δ = T(π) / n

    θ = if isdisc(sphere)
        # --- each face up to half a step before the rim, where the Neumann trace is singular
        face = [(0:(n ÷ 2 - 1)) .* Δ; T(π) / 2 - Δ / 2]
        [face; T(π) .- reverse(face)]
    else
        range(zero(T), T(π); length=n + 1)
    end

    φ = range(zero(T), 2 * T(π); length=round(Int, 360 / resolution) + 1)

    R = frame(sphere)

    # the face of a point on a disc follows from the sign of its coordinate along the axis, which has to survive the
    # rotation into the global frame and back
    δ = sqrt(eps(T)) * equatorialRadius(sphere)

    return [R * surfacePoint(sphere, cos(t), p, δ) for t in θ, p in φ]
end

surfacePoint(sphere::Spheroid, η, φ, δ) = cartesianCoordinates(sphere, SVector(sphere.ξ₀, η, φ))

function surfacePoint(sphere::OblateSpheroid, η, φ, δ)

    point = cartesianCoordinates(sphere, SVector(sphere.ξ₀, η, φ))

    isdisc(sphere) || return point

    return SVector(point[1], point[2], copysign(δ, η))
end


"""
    plotsurface(sphere::Scatterer, values, points; type="abs", scale="linear", normalize=false, clip=1.0, gap=0.02, title="")
    plotsurface(sphere::Scatterer{<:AcousticBoundary}, excitation::AcousticExcitation, Q; total=true, resolution=5, parameter=Parameter(), kwargs...)

Plot a density on the surface of the scatterer as a coloured surface, available after loading PlotlyJS.

The first method plots the `values` given at the `points` of [`surfacePoints`](@ref), which have been computed by
the general API, e.g., `field(sp, ex, PressureTrace(points))`. Optional arguments:
 - `type`: "abs", "real", "imag", or "phase" (in degrees, with a cyclic colour scale)
 - `scale`: "linear" or "log" (in dB, for `type="abs"` only)
 - `normalize`: whether the values are normalized to their largest magnitude
 - `clip`: the quantile of the magnitudes at which the colour scale saturates, e.g., `0.98` for a density that is
   singular at the rim of a disc
 - `gap`: the visual distance of the two faces of a disc, relative to its radius
 - `title`: the title of the colour bar

The second method computes the surface density itself, for `Q` being `PressureTrace`, `PressureNormalGradient`,
or `PressureJump` (a disc only), and plots it: the total density, or the scattered one for `total=false`. The
outward normals of a spheroid are employed for the Neumann trace, and its modal coefficients are computed once.
The keyword arguments `kwargs` are passed on to the first method.

The two faces of a disc are plotted as two surfaces, separated by the `gap` along its axis.
"""
plotsurface(sphere::Scatterer, values, points)
