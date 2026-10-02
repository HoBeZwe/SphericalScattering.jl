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




abstract type Sphere end



"""
    permittivity(sp::Sphere, ex::Excitation, r)

Returns the permittivity of the sphere `sp` at radius `r`.
If this part is PEC, a zero permittivity is returned. The
argument `r` may be a vector of positions.
"""
function permittivity(sp::Sphere, ex::Excitation, r)
    md = medium(sp, ex, r)
    return md.ε
end

"""
    permeability(sp::Sphere, ex::Excitation, r)

Returns the permeability of the sphere `sp` at radius `r`.
If this part is PEC, a zero permeability is returned. The
argument `r` may be a vector of positions.
"""
function permeability(sp::Sphere, ex::Excitation, r)
    md = medium(sp, ex, r)
    return md.μ
end

function permittivity(sp::Sphere, ex::Excitation, pts::AbstractVecOrMat)

    ε = permittivity(sp, ex, norm(first(pts)))
    F = zeros(typeof(ε), size(pts))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(pts)
        F[ind] = permittivity(sp, ex, norm(point))
    end

    return F
end

function permeability(sp::Sphere, ex::Excitation, pts::AbstractVecOrMat)

    ε = permeability(sp, ex, norm(first(pts)))
    F = zeros(typeof(ε), size(pts))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(pts)
        F[ind] = permeability(sp, ex, norm(point))
    end

    return F
end
