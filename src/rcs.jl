
"""
    rcs(sphere::Sphere, excitation::PlaneWave, point_cart; parameter::Parameter=Parameter())

Compute the bistatic radar cross-section (RCS). 
"""
function rcs(sphere::Sphere, excitation::PlaneWave, points_cart; parameter::Parameter=Parameter())

    FF = scatteredfield(sphere, excitation, FarField(points_cart); parameter=parameter)

    return 4π * norm.(FF) .^ 2 / abs2(excitation.amplitude)
end



"""
    rcs(sphere::Sphere, excitation::PlaneWave; parameter::Parameter=Parameter())

Compute the monostatic radar cross-section (RCS): the bistatic RCS solely for the incident direction of the plane wave. 
"""
function rcs(sphere::Sphere, excitation::PlaneWave; parameter::Parameter=Parameter())

    point_cart = -excitation.direction

    return rcs(sphere, excitation, [point_cart]; parameter=parameter)[1]
end



"""
    rcs(sphere::Sphere, excitation::Excitation, point_cart; parameter::Parameter=Parameter())

RCS only defined for plane waves, so far.
"""
function rcs(sphere::Sphere, excitation::Excitation, point_cart; parameter::Parameter=Parameter())

    return error("The (bistatic) RCS is only defined for a plane-wave excitation (so far).")
end



"""
    rcs(sphere::Sphere, excitation::Excitation; parameter::Parameter=Parameter())

RCS only defined for plane waves, so far.
"""
function rcs(sphere::Sphere, excitation::Excitation; parameter::Parameter=Parameter())

    return error("The (monostatic) RCS is only defined for a plane-wave excitation (so far).")
end
