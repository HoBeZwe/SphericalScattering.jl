
# The electromagnetic spheres: a `Sphere` carrying one of the conditions in `boundaries.jl`.
#
# `PECSphere` and `DielectricSphere` are aliases, so that they can be dispatched on as well; the parameters of
# the latter are ordered as those of the former type of the same name. As for the acoustic spheres, the
# docstrings are attached to the constructors rather than to the aliases. The positional constructors are
# kept, as they are in use

const PECSphere = Sphere{PEC}

"""
    PECSphere(
        radius = error("missing argument `radius`"),
    )

Constructor for the PEC sphere, that is, for a [`Sphere`](@ref) with the boundary condition [`PEC`](@ref).

`PECSphere` is an alias of `Sphere{PEC}`, so that it can be dispatched on as well.
"""
PECSphere(; radius=error("missing argument `radius`")) = PECSphere(radius)

PECSphere(radius) = Sphere(radius, PEC())


const DielectricSphere{C,R} = Sphere{Dielectric{C},R}

"""
    DielectricSphere(
        radius      = error("missing argument `radius`"),
        filling     = error("missing argument `filling`")
    )

Constructor for the dielectric sphere, that is, for a [`Sphere`](@ref) with the boundary condition
[`Dielectric`](@ref).

`DielectricSphere` is an alias of `Sphere{<:Dielectric}`, so that it can be dispatched on as well.
"""
DielectricSphere(; radius=error("missing argument `radius`"), filling=error("missing argument `filling`")) =
    DielectricSphere(radius, filling)

DielectricSphere(radius, filling::Medium) = Sphere(radius, Dielectric(filling))



# The sphere with a thin impedance layer and the layered spheres follow the pattern of `PECSphere` and
# `DielectricSphere`. The parameters of the first are ordered as those of the former type of the same name,
# which is not possible for the layered spheres: they store one radius less, the outermost being the radius
# of the sphere

const DielectricSphereThinImpedanceLayer{R,C} = Sphere{ThinImpedanceLayer{R,C},R}

function DielectricSphereThinImpedanceLayer(r::R1, d::R2, thinlayer::Medium{C1}, filling::Medium{C2}) where {R1,R2,C1,C2}

    R = promote_type(R1, R2)
    C = promote_type(C1, C2)

    return Sphere(R(r), ThinImpedanceLayer(R(d), Medium{C}(thinlayer), Medium{C}(filling)))
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

`DielectricSphereThinImpedanceLayer` is an alias of `Sphere{<:ThinImpedanceLayer}`, so that it can be dispatched
on as well, see [`ThinImpedanceLayer`](@ref).
"""
DielectricSphereThinImpedanceLayer(;
    radius=error("missing argument `radius`"),
    thickness=error("missing argument `thickness` of the coating"),
    thinlayer=error("missing argument `thinlayer`"),
    filling=error("missing argument `filling`"),
) = DielectricSphereThinImpedanceLayer(radius, thickness, thinlayer, filling)



const LayeredSphere = Sphere{<:Layered{<:Dielectric}}

"""
    LayeredSphere(
        radii   = error("Missing argument `radii`"),
        filling = error("`missing argument `filling`")
    )

Constructor for the layered dielectric sphere.

`LayeredSphere` is an alias of `Sphere{<:Layered{<:Dielectric}}`, so that it can be dispatched on as well: the
outermost of the `radii` is the radius of the sphere, and the innermost `filling` that of its dielectric core,
see [`Layered`](@ref).
"""
LayeredSphere(; radii=error("Missing argument `radii`"), filling=error("`missing argument `filling`")) = LayeredSphere(radii, filling)

function LayeredSphere(radii::SVector, filling::SVector)

    if sort(radii) != radii
        error("Radii are not ordered ascendingly.")
    end

    if length(radii) != length(filling)
        error("Number of fillings does not match number of radii.")
    end

    return Sphere(last(radii), Layered(pop(radii), popfirst(filling), Dielectric(first(filling))))
end



const LayeredSpherePEC = Sphere{<:Layered{PEC}}

"""
    LayeredSpherePEC(
        radii   = error("Missing argument `radii`"),
        filling = error("Missing argument `filling`")
    )

Constructor for the layered dielectric sphere with a PEC core.

`LayeredSpherePEC` is an alias of `Sphere{<:Layered{PEC}}`, so that it can be dispatched on as well: the
outermost of the `radii` is the radius of the sphere, and the innermost one that of its PEC core, see
[`Layered`](@ref).
"""
LayeredSpherePEC(; radii=error("Missing argument `radii`"), filling=error("Missing argument `filling`")) =
    LayeredSpherePEC(radii, filling)

function LayeredSpherePEC(radii::SVector, filling::SVector)

    if sort(radii) != radii
        error("Radii are not ordered ascendingly.")
    end

    if length(radii) != length(filling) + 1
        error("Number of fillings does not match number of radii.")
    end

    return Sphere(last(radii), Layered(pop(radii), filling, PEC()))
end



"""
    layerRadii(sp::Sphere{<:Layered})

Returns the radii of all interfaces of a layered sphere ascendingly, the outermost being the radius of the sphere.
"""
layerRadii(sp::Sphere{<:Layered}) = push(sp.boundary.radii, sp.radius)

"""
    layerFillings(sp::Sphere{<:Layered})

Returns the fillings of a layered sphere from the inside out: those of the shells, preceded by that of the core
if it is a dielectric one.
"""
layerFillings(sp::Sphere{<:Layered{<:Dielectric}}) = pushfirst(sp.boundary.filling, sp.boundary.core.filling)

layerFillings(sp::Sphere{<:Layered{PEC}}) = sp.boundary.filling



"""
    numlayers(sp::Scatterer{<:ElectromagneticBoundary})

Returns the number of layers.
"""
function numlayers(sp::Scatterer{<:ElectromagneticBoundary})
    return 2
end

function numlayers(sp::Sphere{<:Layered})
    return length(layerRadii(sp)) + 1
end



"""
    layer(sp::Scatterer{<:ElectromagneticBoundary}, r)

Returns the index of the layer, `r` is located, where `1` denotes the inner most layer. 
"""
function layer(sp::Scatterer{<:ElectromagneticBoundary}, r)
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
    layer(sp::Sphere{<:Layered}, r)

Returns the index of the layer, `r` is located, where `1` denotes the inner most layer.
"""
function layer(sp::Sphere{<:Layered}, r)
    # Using Jin's numbering from the multi-layered cartesian
    # in anticipation. For PEC and dielectric sphere;
    # 1 = interior, 2 = exterior

    r < 0.0 && error("The radius must be a positive number.")

    N = numlayers(sp)
    radii = layerRadii(sp)

    for i in 1:(N - 1)
        if r < radii[i]
            return i # Convention: Boundary belongs to outer layer
        end
    end

    return N
end


"""
wavenumber(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)

Returns the wavenumber at radius `r` in the sphere `sp`.
If this part is PEC, k=0 is returned.
"""
function wavenumber(sp::Sphere{PEC}, ex::Excitation, r)
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

function wavenumber(sp::Sphere{<:Dielectric}, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.boundary.filling.ε
        μ = sp.boundary.filling.μ
    end

    c = 1 / sqrt(ε * μ)
    k = 2π * ex.frequency / c

    return k
end

"""
    impedance(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)

Returns the wavenumber at radius `r` in the sphere `sp`.
If this part is PEC, a zero medium is returned.
"""
function impedance(sp::Sphere{<:Dielectric}, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        ε = sp.boundary.filling.ε
        μ = sp.boundary.filling.μ
    end

    return sqrt(μ / ε)
end

"""
    impedance(sp::Scatterer{<:ElectromagneticBoundary}, r)

Returns the impedance of the sphere `sp` at radius `r`.
If this part is PEC, a zero medium is returned.
"""
function impedance(sp::Sphere{PEC}, ex::Excitation, r)
    if layer(sp, r) == 2
        ε = ex.embedding.ε
        μ = ex.embedding.μ
    else
        return promote_type(typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
    end

    return sqrt(μ / ε)
end

function medium(sp::Sphere{<:Layered{<:Dielectric}}, ex::Excitation, r)
    N = numlayers(sp) # Number of interior layers

    if layer(sp, r) == N # Outer layer has largest index
        return ex.embedding
    else
        return layerFillings(sp)[layer(sp, r)]
    end
end



function medium(sp::Sphere{<:Layered{PEC}}, ex::Excitation, r)
    N = numlayers(sp) # Number of interior layers

    fillings = layerFillings(sp)

    if layer(sp, r) == N # Outer layer has largest index
        return ex.embedding
    elseif layer(sp, r) == 1
        Z = promote_type(typeof(fillings[end].ε), typeof(fillings[end].μ), typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
        return Medium(Z, Z)
    else
        return fillings[layer(sp, r) - 1]
    end
end



function medium(sp::Sphere{<:Union{Dielectric,ThinImpedanceLayer}}, ex::Excitation, r)
    if layer(sp, r) == 2
        return ex.embedding
    else
        return sp.boundary.filling
    end
end



"""
    medium(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)

Returns the medium of the sphere `sp` at radius `r`.
If this part is PEC, a zero medium is returned. The
argument `r` may be a vector of positions.
"""
function medium(sp::Sphere{PEC}, ex::Excitation, r)
    if layer(sp, r) == 2
        return ex.embedding
    else
        Z = promote_type(typeof(ex.embedding.ε), typeof(ex.embedding.μ))(0.0)
        return Medium(Z, Z)
    end
end



"""
    permittivity(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)

Returns the permittivity of the sphere `sp` at radius `r`.
If this part is PEC, a zero permittivity is returned. The
argument `r` may be a vector of positions.
"""
function permittivity(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)
    md = medium(sp, ex, r)
    return md.ε
end

"""
    permeability(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)

Returns the permeability of the sphere `sp` at radius `r`.
If this part is PEC, a zero permeability is returned. The
argument `r` may be a vector of positions.
"""
function permeability(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, r)
    md = medium(sp, ex, r)
    return md.μ
end

function permittivity(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, pts::AbstractVecOrMat)

    ε = permittivity(sp, ex, norm(first(pts)))
    F = zeros(typeof(ε), size(pts))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(pts)
        F[ind] = permittivity(sp, ex, norm(point))
    end

    return F
end

function permeability(sp::Scatterer{<:ElectromagneticBoundary}, ex::Excitation, pts::AbstractVecOrMat)

    ε = permeability(sp, ex, norm(first(pts)))
    F = zeros(typeof(ε), size(pts))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(pts)
        F[ind] = permeability(sp, ex, norm(point))
    end

    return F
end
