

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

abstract type Excitation end

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
