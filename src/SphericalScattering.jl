module SphericalScattering


const μ0 = 4pi * 1e-7        # default permeability
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
export ElectricRingCurrent, MagneticRingCurrent
export FarField, ElectricField, MagneticField
export Medium, Parameter

# functions
export electricRingCurrent, magneticRingCurrent
export HertzianDipole, FitzgeraldDipole
export planeWave
export SphericalModeTE, SphericalModeTM
export PECSphere
export field, scatteredfield



# -------- included files
include("coordinateTransforms.jl")
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
end
