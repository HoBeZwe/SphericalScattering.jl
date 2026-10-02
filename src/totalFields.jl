
inside(sphere::Sphere) = 0.0
inside(sphere::PECSphere) = sphere.radius - 1e-15
inside(sphere::Union{HardSphere,SoftSphere}) = sphere.radius - 1e-15


"""
    field(sphere::Sphere, excitation::Excitation, quantity::Field; parameter::Parameter=Parameter())

Compute the total field in the presence of a sphere for a given excitation.
"""
function field(sphere::Sphere, excitation::Excitation, quantity::Field; parameter::Parameter=Parameter())

    # incident and scattered field
    F = field(excitation, quantity; parameter=parameter, zeroRadius=inside(sphere))
    F .+= scatteredfield(sphere, excitation, quantity; parameter=parameter)

    return F
end



"""
    field(sphere::Sphere, excitation::Excitation, quantity::Trace; parameter::Parameter=Parameter())

Compute the total trace on the surface of a sphere for a given excitation.

In contrast to the fields, the traces are not set to zero anywhere: they are defined on the surface of the sphere,
where the locations are assumed to lie.
"""
function field(sphere::Sphere, excitation::Excitation, quantity::Trace; parameter::Parameter=Parameter())

    # incident and scattered trace
    F = field(excitation, quantity; parameter=parameter)
    F .+= scatteredfield(sphere, excitation, quantity; parameter=parameter)

    return F
end



"""
    field(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Descriptive error for the total far-field in the presence of a sphere for an incident plane wave.
"""
function field(sphere::Sphere, excitation::PlaneWave, quantity::FarField; parameter::Parameter=Parameter())

    return error("The total far-field for a plane-wave excitation is not defined")
end


"""
    field(sphere::Sphere, excitation::AcousticPlaneWave, quantity::Field; parameter::Parameter=Parameter())

Descriptive error for the total far-field in the presence of a sphere for an incident acoustic plane wave.
"""
function field(sphere::Sphere, excitation::AcousticPlaneWave, quantity::FarField; parameter::Parameter=Parameter())

    return error("The total far-field for a plane-wave excitation is not defined")
end


"""
    field(sphere::Sphere, excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

Descriptive error for the total far-field in the presence of a sphere for an incident spherical mode.
"""
function field(sphere::Sphere, excitation::SphericalMode, quantity::FarField; parameter::Parameter=Parameter())

    return error("The total far-field for a spherical mode excitation is not defined")
end
