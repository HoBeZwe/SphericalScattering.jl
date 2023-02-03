
"""
    field(sphere::Sphere, excitation::Excitation, quantity::Field; parameter::Parameter=Parameter())

Compute the total field in the presence of a sphere for a given excitation.
"""
function field(sphere::Sphere, excitation::Excitation, quantity::Field; parameter::Parameter=Parameter())

    # incident and scattered field
    Finc = field(excitation, quantity; parameter=parameter)
    Fsca = scatteredfield(sphere, excitation, quantity; parameter=parameter)

    return Finc .+ Fsca
end



"""
    field(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Descriptive error for the total far-field in the presence of a sphere for an incident plane wave.
"""
function field(sphere::Sphere, excitation::PlaneWave, quantity::FarField; parameter::Parameter=Parameter())

    return error("The total far-field for a plane-wave excitation is not defined")
end


"""
    field(sphere::Sphere, excitation::SphericalMode, quantity::Field; parameter::Parameter=Parameter())

Descriptive error for the total far-field in the presence of a sphere for an incident spherical mode.
"""
function field(sphere::Sphere, excitation::SphericalMode, quantity::FarField; parameter::Parameter=Parameter())

    return error("The total far-field for a spherical mode excitation is not defined")
end
