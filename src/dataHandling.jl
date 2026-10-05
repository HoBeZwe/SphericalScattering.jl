

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
    PressureJump(locations)

The jump ``[p] = p|_+ - p|_-`` of the pressure across an open surface at the given locations.

The faces are distinguished by the outward normal: ``+`` denotes the one whose normal is the axis of the
scatterer. The jump is the natural unknown of a boundary element formulation on an open surface, where a
one-sided trace is not determined by the geometry alone.
"""
struct PressureJump <: Trace
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

# the excitations are split by physics, mirroring the boundary conditions of the scatterers, so that a
# scatterer and an excitation of different physics can be rejected by dispatch
abstract type ElectromagneticExcitation <: Excitation end

# the acoustic excitations share the series for the scattered field, which differs only in the axis of
# rotational symmetry and in the coefficients of the incident expansion
abstract type AcousticExcitation <: Excitation end

wavenumber(ex::Excitation) = 2π * ex.frequency * sqrt(ex.embedding.ε * ex.embedding.μ)

#abstract type Parameter end

struct Parameter
    nmax::Int
    relativeAccuracy::AbstractFloat
end

"""
    Parameter(nmax = -1, relativeAccuracy = 1e-12)

Settings controlling the truncation of the series which are evaluated.

The `relativeAccuracy` is the relative contribution below which a term is considered negligible: series which
converge on their own are terminated once a term contributes less than that, and truncations which are
determined in advance are chosen such that the neglected terms stay below it.

A non-negative `nmax` fixes the truncation instead, overriding the automatic choice; the default of `-1` leaves
it to the implementation.
"""
Parameter() = Parameter(-1, 1e-12)

# global setting for the style of the progress bar
function progress(numIter::Int)
    return Progress(numIter; barglyphs=BarGlyphs("[=> ]"), color=:white)
end
