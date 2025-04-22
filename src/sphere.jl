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

struct DielectricSphere{C,R} <: Sphere
    radius::R
    filling::Medium{C}
end


"""
    DielectricSphere(
        radius      = error("missing argument `radius`"), 
        filling     = error("missing argument `filling`")
    )

Constructor for the dielectric sphere.
"""
DielectricSphere(; radius=error("missing argument `radius`"), filling=error("missing argument `filling`")) = DielectricSphere(
    radius, filling
)



struct DielectricSphereThinImpedanceLayer{R,C} <: Sphere
    radius::R
    thickness::R
    thinlayer::Medium{C}
    filling::Medium{C}
end

function DielectricSphereThinImpedanceLayer(r::R1, d::R2, thinlayer::Medium{C1}, filling::Medium{C2}) where {R1,R2,C1,C2}

    R = promote_type(R1, R2)
    C = promote_type(C1, C2)

    DielectricSphereThinImpedanceLayer(R(r), R(d), Medium{C}(thinlayer), Medium{C}(filling))
end

"""
    DielectricSphereThinImpedanceLayer(
        radius    = error("missing argument `radius`"),
        thickness = error("missing argument `thickness` of the coating"),
        thinlayer = error("missing argument `thinlayer`"),
        filling   = error("missing argument `filling`")
    )

Constructor for the dielectric sphere with a thin impedance layer.
For this model, it is assumed that the displacement field is only radial
direction in the layer, which requires a small thickness and low conductivity.
For details, see for example T. B. Jones, Ed., “Models for layered spherical particles,”
in Electromechanics of Particles, Cambridge: Cambridge University Press, 1995, 
pp. 227–235. doi: 10.1017/CBO9780511574498.012.
"""
DielectricSphereThinImpedanceLayer(; radius=error("missing argument `radius`"), thickness=error("missing argument `thickness` of the coating"), thinlayer=error("missing argument `thinlayer`"), filling=error("missing argument `filling`")) = DielectricSphereThinImpedanceLayer(
    radius, thickness, thinlayer, filling
)



struct PECSphere{R} <: Sphere
    radius::R
end

"""
    PECSphere( 
        radius = error("missing argument `radius`"), 
    )

Constructor for the PEC sphere.
"""
PECSphere(; radius=error("missing argument `radius`")) = PECSphere(radius)



struct LayeredSphere{N,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{N,Medium{C}}
end

"""
    LayeredSphere( 
        radii   = error("Missing argument `radii`"), 
        filling = error("`missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
function LayeredSphere(; radii=error("Missing argument `radii`"), filling=error("`missing argument `filling`"))

    if sort(radii) != radii
        error("Radii are not ordered ascendingly.")
    end

    if length(radii) != length(filling)
        error("Number of fillings does not match number of radii.")
    end

    LayeredSphere(radii, filling)
end



struct LayeredSpherePEC{N,D,R,C} <: Sphere
    radii::SVector{N,R}
    filling::SVector{D,Medium{C}}
end

"""
    LayeredSpherePEC( 
        radii   = error("Missing argument `radii`"), 
        filling = error("Missing argument `filling`")
    )

Constructor for the layered dielectric sphere.
"""
function LayeredSpherePEC(; radii=error("Missing argument `radii`"), filling=error("Missing argument `filling`"))

    if sort(radii) != radii
        error("Radii are not ordered ascendingly.")
    end

    if length(radii) != length(filling) + 1
        error("Number of fillings does not match number of radii.")
    end

    LayeredSpherePEC(radii, filling)
end


"""
    numlayers(sp::Sphere)

Returns the number of layers.
"""
function numlayers(sp::Sphere)
    return 2
end

function numlayers(sp::Union{LayeredSphere,LayeredSpherePEC})
    return length(sp.radii) + 1
end



"""
    layer(sp::Sphere, r)

Returns the index of the layer, `r` is located, where `1` denotes the inner most layer. 
"""
function layer(sp::Sphere, r)
    # Using Jin's numbering from the multi-layered cartesian
    # in anticipation. For PEC and dielectric sphere;
    # 1 = interior, 2 = exterior

    if r >= sp.radius
        return 2
    else
        return 1
    end
end



"""
    layer(sp::LayeredSphere, r)

Returns the index of the layer, `r` is located, where `1` denotes the inner most layer. 
"""
function layer(sp::Union{LayeredSphere,LayeredSpherePEC}, r)
    # Using Jin's numbering from the multi-layered cartesian
    # in anticipation. For PEC and dielectric sphere;
    # 1 = interior, 2 = exterior

    r < 0.0 && error("The radius must be a positive number.")

    N = numlayers(sp)

    for i in 1:(N - 1)
        if r < sp.radii[i]
            return i # Convention: Boundary belongs to outer layer
        end
    end

    return N
end


"""
wavenumber(sp::Sphere, ex::Excitation, r)

Returns the wavenumber at radius `r` in the sphere `sp`.
If this part is PEC, k=0 is rweturned.
"""
function wavenumber(sp::PECSphere, ex::Excitation, r)
    ε = ex.embedding.ε
    μ = ex.embedding.μ

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    if layer(sp, r) == 2
        return k
    else
        return typeof(k)(0.0)
    end
end

function wavenumber(sp::DielectricSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    return k
end

"""
    impedance(sp::Sphere, ex::Excitation, r)

Returns the wavenumber at radius `r` in the sphere `sp`.
If this part is PEC, a zero medium is returned.
"""
function impedance(sp::DielectricSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.filling.ε
        μ = sp.filling.μ
    end

    return sqrt(μ / ε)
end

"""
    impedance(sp::Sphere, r)

Returns the impedance of the sphere `sp` at radius `r`.
If this part is PEC, a zero medium is returned.
"""
function impedance(sp::PECSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        return promote_type(typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
    end

    return sqrt(μ / ε)
end

function medium(sp::LayeredSphere, ex::Excitation, r)
    N = numlayers(sp) # Number of interior layers

    if layer(sp, r) == N # Outer layer has largest index
        return ex.embedding
    else
        return sp.filling[layer(sp, r)]
    end
end



function medium(sp::LayeredSpherePEC, ex::Excitation, r)
    N = numlayers(sp) # Number of interior layers

    if layer(sp, r) == N # Outer layer has largest index
        return ex.embedding
    elseif layer(sp, r) == 1
        Z = promote_type(typeof(sp.filling[end].ε), typeof(sp.filling[end].μ), typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
        return Medium(Z, Z)
    else
        return sp.filling[layer(sp, r) - 1]
    end
end



function medium(sp::Union{DielectricSphere,DielectricSphereThinImpedanceLayer}, ex::Excitation, r)
    if layer(sp, r) == 2
        return ex.embedding
    else
        return sp.filling
    end
end



"""
    medium(sp::Sphere, ex::Excitation, r)

Returns the medium of the sphere `sp` at radius `r`.
If this part is PEC, a zero medium is returned. The
argument `r` may be a vector of positions.
"""
function medium(sp::PECSphere, ex::Excitation, r)
    if layer(sp, r) == 2
        return ex.embedding
    else
        Z = promote_type(typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
        return Medium(Z, Z)
    end
end

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
