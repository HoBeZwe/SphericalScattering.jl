
# the quantities obtained from the series for the scattered pressure itself; the Neumann trace
# requires the derivatives of that series and is handled separately
const AcousticQuantity = Union{Pressure,FarField,PressureTrace}


"""
    checkExcitation(scatterer::Scatterer{<:AcousticBoundary}, excitation::AcousticExcitation)

Ensure that the excitation is compatible with the scatterer; nothing has to be checked by default.

Whether the scatterer is an acoustic one at all is decided by dispatch instead, see [`Scatterer`](@ref).
"""
checkExcitation(scatterer::Scatterer{<:AcousticBoundary}, excitation::AcousticExcitation) = nothing
