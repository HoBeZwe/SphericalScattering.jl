

abstract type Field end

struct FarField <: Field
    locations#::Vector
end

struct ElectricField <: Field
    locations#::Vector
end

struct DisplacementField <: Field
    locations#::Vector
end

struct MagneticField <: Field
    locations#::Vector
end

struct ScalarPotential <: Field
    locations
end

struct ScalarPotentialJump <: Field
    locations
end

struct Pressure <: Field
    locations
end



abstract type Trace end

struct PressureTrace <: Trace
    locations
end

"""
    PressureNormalGradient(locations)
    PressureNormalGradient(locations, normals)

The normal gradient ``\\hat{n} ⋅ ∇p`` of the pressure at the given locations.

If no `normals` are provided, the outward normal ``\\hat{n} = \\hat{r}`` of a sphere centered at the origin is
employed at every location. Otherwise one normal vector per location is expected, which is useful for the
locations of a faceted surface mesh, whose normals do not coincide with ``\\hat{r}``. Either way, the normals are
normalized and stored for every location.
"""
struct PressureNormalGradient <: Trace
    locations
    normals

    # inner constructor: default to the outward normal n̂ = r̂ and normalize
    function PressureNormalGradient(locations, normals=nothing)

        isnothing(normals) && return new(locations, normalize.(locations))

        length(normals) == length(locations) || error("The number of provided normal vectors does not match the number of locations.")

        return new(locations, normalize.(normals))
    end
end



abstract type Excitation end

# the acoustic excitations share the series for the scattered field, which differs only in the axis of
# rotational symmetry and in the coefficients of the incident expansion
abstract type AcousticExcitation <: Excitation end

wavenumber(ex::Excitation) = 2π * ex.frequency * sqrt(ex.embedding.ε * ex.embedding.μ)

#abstract type Parameter end

struct Parameter
    nmax::Int
    relativeAccuracy::AbstractFloat
end

Parameter() = Parameter(-1, 1e-12)

# global setting for the style of the progress bar
function progress(numIter::Int)
    return Progress(numIter; barglyphs=BarGlyphs("[=> ]"), color=:white)
end
