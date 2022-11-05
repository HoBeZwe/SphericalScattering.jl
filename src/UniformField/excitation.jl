struct UniformField{R,C,T} <: Excitation
    embedding::Medium{C}
    amplitude::T
    direction::SVector{3,R}
end

UniformField(; embedding=Medium(ε0, μ0), amplitude=1.0, direction=SVector{3,Float64}(1.0, 0.0, 0.0)) =
    UniformField(embedding, amplitude, direction)
