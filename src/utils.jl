
"""
    patternPoints(;r=1.0, resolution=5)

Convenience function returning points (in Cartesian and spherical coordinates) on a spherical grid at a distance 'r' and a resolution in degrees. 
"""
function sphericalGridPoints(;r=1.0, resolution=5)

    ϑ = range(0.0; stop=π, length=round(Int, 180 / resolution))  # default: 5° steps
    ϕ = range(0.0; stop=2π, length=round(Int, 360 / resolution)) # default: 5° steps 

    points_cart = [SVector(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]
    points_sph  = [SVector(r, θ, φ) for θ in ϑ, φ in ϕ]

    return points_cart, points_sph
end