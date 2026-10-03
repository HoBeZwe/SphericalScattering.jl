
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


"""
    checkScatterer(scatterer::Sphere)

Ensure that an acoustic solution is implemented for the scatterer.

The check belongs before the loop over the locations: an error thrown inside the parallel loop is wrapped in a
`TaskFailedException` as soon as more than one thread is available, which would make the failure depend on the
number of threads.
"""
checkScatterer(scatterer::AcousticScatterer) = nothing

checkScatterer(scatterer::Sphere) = error("Acoustic scattering is only implemented for sound-hard and sound-soft spheres (so far).")
