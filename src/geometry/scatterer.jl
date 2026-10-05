
# The scatterers are the product of a geometry and a condition on its surface. The geometry is the type,
# the condition its parameter, so that methods can dispatch on either alone or on both: on `Sphere` for a
# geometric query, on `Sphere{SoundHard}` for a solution, and on `Scatterer{<:AcousticBoundary}` for
# what all scatterers of one physics share.


"""
    Boundary

Supertype of the conditions on the surface of a scatterer.

The conditions are orthogonal to the geometry of the scatterer and are, therefore, carried as its type parameter.
They are split by physics into [`AcousticBoundary`](@ref) and [`ElectromagneticBoundary`](@ref).
"""
abstract type Boundary end

"""
    AcousticBoundary

Supertype of the acoustic boundary conditions, see [`SoundHard`](@ref) and [`SoundSoft`](@ref).
"""
abstract type AcousticBoundary <: Boundary end

"""
    ElectromagneticBoundary

Supertype of the electromagnetic boundary conditions, see [`PEC`](@ref), [`Dielectric`](@ref), [`Layered`](@ref)
and [`ThinImpedanceLayer`](@ref).
"""
abstract type ElectromagneticBoundary <: Boundary end

# the conditions themselves are defined with their physics, in `acoustics/boundaries.jl` and
# `electromagnetics/boundaries.jl`



"""
    Scatterer{BC}

Supertype of the scatterers, `BC` being the condition on their surface, see [`Boundary`](@ref).

The parameter allows to address all scatterers of one physics at once, e.g., `Scatterer{<:AcousticBoundary}`.
Solutions are, however, dispatched on a concrete geometry and condition, such as `Sphere{SoundHard}`; the
abstract form is meant for preconditions and descriptive errors, lest the methods become ambiguous.
"""
abstract type Scatterer{BC<:Boundary} end



"""
    Sphere{BC,R} <: Scatterer{BC}

Sphere centered in the origin, on whose surface the condition `boundary` of type `BC` holds.

The condition is stored as a value besides being the type parameter: a parameter-free condition, such as
[`SoundHard`](@ref), occupies no memory, whereas a condition requiring data can carry it.

Not to be confused with a [`Spheroid`](@ref), which it is not a special case of: the spheroidal coordinates
degenerate for a sphere.
"""
struct Sphere{BC<:Boundary,R} <: Scatterer{BC}
    radius::R
    boundary::BC
end

"""
    Sphere(
        radius   = error("missing argument `radius`"),
        boundary = error("missing argument `boundary`")
    )

Constructor for a sphere of the given `radius`, on whose surface the condition `boundary` holds.
"""
Sphere(; radius=error("missing argument `radius`"), boundary=error("missing argument `boundary`")) = Sphere(radius, boundary)

"""
    Sphere{BC}(
        radius = error("missing argument `radius`")
    )

Constructor for a sphere with a parameter-free boundary condition, e.g., `Sphere{SoundHard}(; radius=1.0)`.

A condition requiring data has to be passed as a value instead, see [`Sphere`](@ref).
"""
function Sphere{BC}(; radius=error("missing argument `radius`")) where {BC<:Boundary}

    Base.issingletontype(BC) || error(
        "`Sphere{BC}` requires a parameter-free boundary condition, which $BC is not: pass the condition as a value, `Sphere(; radius, boundary)`.",
    )

    return Sphere(radius, BC())
end



"""
    isinside(scatterer::Scatterer, point)

Returns whether the point lies inside the scatterer.

Which geometry the interior has is for the scatterer to answer: the radius of a sphere cannot express the
interior of a spheroid, nor the other way around. The methods are therefore found with the respective types.
"""
function isinside end

"""
    isinside(scatterer::Sphere, point)

Returns whether the point lies inside the sphere, see [`isinside`](@ref).

Being a question of geometry alone, it holds for every boundary condition.
"""
isinside(scatterer::Sphere, point) = norm(point) < scatterer.radius



# --- a scatterer and an excitation of different physics. These are checked on the level of the locations,
#     by dispatch, so that the error is thrown before any loop over them: inside a parallel loop it would be
#     wrapped in a `TaskFailedException` as soon as more than one thread is available

"""
    scatteredfield(sphere::Scatterer{<:AcousticBoundary}, excitation::ElectromagneticExcitation, quantity)

Descriptive error for an acoustic scatterer excited by an electromagnetic field.
"""
function scatteredfield(sphere::Scatterer{<:AcousticBoundary}, excitation::ElectromagneticExcitation, quantity; kwargs...)

    return error(
        "An electromagnetic excitation requires a scatterer with an electromagnetic boundary condition, such as a `PECSphere`."
    )
end

"""
    scatteredfield(sphere::Scatterer{<:ElectromagneticBoundary}, excitation::AcousticExcitation, quantity)

Descriptive error for an electromagnetic scatterer excited by an acoustic field.
"""
function scatteredfield(sphere::Scatterer{<:ElectromagneticBoundary}, excitation::AcousticExcitation, quantity; kwargs...)

    return error(
        "An acoustic excitation requires a scatterer with an acoustic boundary condition, such as a `HardSphere`, a `SoftSphere` or a `Spheroid`.",
    )
end
