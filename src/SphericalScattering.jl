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
using StaticArrays
using OhMyThreads
using ProgressMeter



# -------- exportet parts
# types
export Excitation
export PlaneWave
export UniformField
export ElectricRingCurrent, MagneticRingCurrent
export FarField, ElectricField, MagneticField
export DisplacementField
export ScalarPotential, ScalarPotentialJump
export Medium, Parameter
export μ0, ε0

export Acoustic
export AcousticPlaneWave, AcousticMonopole
export Pressure
export PressureTrace, PressureNormalGradient

# functions
export electricRingCurrent, magneticRingCurrent
export HertzianDipole, FitzgeraldDipole
export planeWave
export SphericalMode, SphericalModeTE, SphericalModeTM
export PECSphere, DielectricSphere, LayeredSphere, LayeredSpherePEC
export DielectricSphereThinImpedanceLayer
export HardSphere, SoftSphere
export field, scatteredfield
export rcs
export sphericalGridPoints, phiCutPoints, thetaCutPoints
export numlayers, layer
export permittivity, permeability, medium


# -------- extensions
function plotff end
function plotnf end
function plotffcut end
function plotnfcut end

export plotff, plotnf, plotffcut, plotnfcut


# -------- included files
include("dataHandling.jl")
include("sphere.jl")

include("electromagnetics/sphere.jl")

include("electromagnetics/ringCurrent/excitation.jl")
include("electromagnetics/ringCurrent/incident.jl")
include("electromagnetics/ringCurrent/scattered.jl")

include("electromagnetics/dipoles/excitation.jl")
include("electromagnetics/dipoles/incident.jl")
include("electromagnetics/dipoles/scattered.jl")

include("electromagnetics/planeWave/excitation.jl")
include("electromagnetics/planeWave/incident.jl")
include("electromagnetics/planeWave/scattered.jl")

include("electromagnetics/sphericalModes/excitation.jl")
include("electromagnetics/sphericalModes/incident.jl")
include("electromagnetics/sphericalModes/scattered.jl")

include("electromagnetics/UniformField/excitation.jl")
include("electromagnetics/UniformField/incident.jl")
include("electromagnetics/UniformField/scattered.jl")


include("acoustics/sphere.jl")
include("acoustics/scattered.jl") # the series shared by all acoustic excitations

include("acoustics/planeWave/excitation.jl")
include("acoustics/planeWave/incident.jl")
include("acoustics/planeWave/scattered.jl")

include("acoustics/monopole/excitation.jl")
include("acoustics/monopole/incident.jl")
include("acoustics/monopole/scattered.jl")

include("acoustics/Acoustic.jl") # the single `Acoustic` submodule holding the user-facing constructors

include("totalFields.jl")
include("coordinateTransforms.jl")
include("utils.jl")
include("rcs.jl")

if !isdefined(Base, :get_extension)
    include("../ext/SphericalScatteringExt.jl") # for backwards compatibility with julia versions below 1.9
end

end
