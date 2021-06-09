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
export ElectricRingCurrent, MagneticRingCurrent
export FarField, ElectricField, MagneticField
export Medium, Parameter

# functions
export electricRingCurrent
export PECSphere
export field, scatteredfield



# -------- included files
include("coordinateTransforms.jl")
include("dataHandling.jl")
include("sphere.jl")

include("ringCurrent/excitation.jl")
include("ringCurrent/incident.jl")
include("ringCurrent/scattered.jl")

end
