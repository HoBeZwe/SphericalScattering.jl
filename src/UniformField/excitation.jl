
struct UniformField{C,T,R} <: Excitation
    embedding::Medium{C}
    amplitude::T
    direction::SVector{3,R}

    # inner constructor: normalize direction
    function UniformField(embedding::Medium{C}, amplitude::T, direction::SVector{3,R}) where {C,R,T}
        dir_normalized = normalize(direction)
        new{C,R,T}(embedding, amplitude, dir_normalized)
    end
end


"""
    ex = UniformField(; 
        embedding=Medium(ε0, μ0), 
        amplitude=1.0, 
        direction=SVector{3,Float64}(1.0, 0.0, 0.0)
    )
"""
UniformField(; embedding=Medium(ε0, μ0), amplitude=1.0, direction=SVector{3,Float64}(1.0, 0.0, 0.0)) =
    UniformField(embedding, amplitude, direction)
