
struct AcousticMonopole{T,R,C} <: Excitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
end


# module Acoustic

# using ..SphericalScattering
# using StaticArrays

# """
#     monopole(;
#         position=error("missing argument `position`"),
#         frequency=error("missing argument `frequency`"),
#         embedding=SphericalScattering.Medium(ε0, μ0),    
#         amplitude=1.0
#     ) 
# """
# monopole(;
#     position=error("missing argument `position`"),
#     frequency=error("missing argument `frequency`"),
#     embedding=SphericalScattering.Medium(ε0, μ0),
#     amplitude=1.0
# ) = SphericalScattering.AcousticMonopole(embedding, frequency, amplitude, position)

# end