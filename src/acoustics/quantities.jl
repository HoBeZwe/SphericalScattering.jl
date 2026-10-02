
# the quantities obtained from the series for the scattered pressure itself; the Neumann trace
# requires the derivatives of that series and is handled separately
const AcousticQuantity = Union{Pressure,FarField,PressureTrace}

# the scatterers for which an acoustic solution is implemented. Scatterers of the electromagnetic part
# are subtypes of `Sphere` as well, so the preconditions have to be restricted to these, lest they
# shadow the descriptive error of the unsupported ones
const AcousticScatterer = Union{HardSphere,SoftSphere,Spheroid}


"""
    isinside(scatterer::Sphere, point)

Returns whether the point lies inside the scatterer.

Which geometry the interior has is for the scatterer to answer: the radius of a sphere cannot express the
interior of a spheroid, nor the other way around. The methods are therefore found with the respective types.
"""
function isinside end


"""
    checkExcitation(scatterer::Sphere, excitation::AcousticExcitation)

Ensure that the excitation is compatible with the scatterer; nothing has to be checked by default.
"""
checkExcitation(scatterer::Sphere, excitation::AcousticExcitation) = nothing
