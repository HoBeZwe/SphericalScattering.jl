"""
    Medium(ε, μ)

Homogeneous, isotropic background medium.
"""
struct Medium{C}
    ε::C
    μ::C
end

function Medium(ε::T1, μ::T2) where {T1,T2}
    T = promote_type(T1, T2)
    return Medium(T(ε), T(μ))
end

function Medium{T}(md) where {T}
    return Medium(T(md.ε), T(md.μ))
end
