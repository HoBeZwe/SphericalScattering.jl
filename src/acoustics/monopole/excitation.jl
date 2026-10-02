
struct AcousticMonopole{T,R,C} <: AcousticExcitation
    embedding::Medium{C}
    frequency::R
    amplitude::T
    position::SVector{3,R}
end
