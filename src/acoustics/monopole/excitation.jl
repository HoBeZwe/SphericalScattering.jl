
struct AcousticMonopole{T,R,C} <: Excitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
end
