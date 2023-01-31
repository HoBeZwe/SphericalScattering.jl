module SphericalScattering

"""
    μ0 = 4pi * 1e-7 

Free space permeability.
"""
const μ0 = 4pi * 1e-7        # default permeability

"""
    ε0 = 8.8541878176e-12 

Free space permittivity.
"""
const ε0 = 8.8541878176e-12  # default permittivity



# -------- used packages
using SpecialFunctions, LegendrePolynomials
using LinearAlgebra
#using Combinatorics
using StaticArrays
using Requires



# -------- exportet parts
# types
export Excitation
export PlaneWave
export UniformField
export ElectricRingCurrent, MagneticRingCurrent
export FarField, ElectricField, MagneticField
export ScalarPotential
export Medium, Parameter
export μ0, ε0

# functions
export electricRingCurrent, magneticRingCurrent
export HertzianDipole, FitzgeraldDipole
export planeWave
export SphericalMode, SphericalModeTE, SphericalModeTM
export PECSphere, DielectricSphere, LayeredSphere, LayeredSpherePEC
export field, scatteredfield



# -------- included files
include("dataHandling.jl")
include("sphere.jl")

include("ringCurrent/excitation.jl")
include("ringCurrent/incident.jl")
include("ringCurrent/scattered.jl")

include("dipoles/excitation.jl")
include("dipoles/incident.jl")
include("dipoles/scattered.jl")

include("planeWave/excitation.jl")
include("planeWave/incident.jl")
include("planeWave/scattered.jl")

include("sphericalModes/excitation.jl")
include("sphericalModes/incident.jl")
include("sphericalModes/scattered.jl")

include("UniformField/excitation.jl")
include("UniformField/incident.jl")
include("UniformField/scattered.jl")

include("totalFields.jl")
include("coordinateTransforms.jl")
end
