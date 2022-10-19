
"""
field(excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field of a uniform field.
"""
function field(excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

F = zeros(fieldType(quantity), length(quantity.locations))

# --- compute field in Cartesian representation
for (ind, point) in enumerate(quantity.locations)
    F[ind] = field(excitation, point, quantity, parameter=parameter)
end

return F
end



"""
field(excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the vector field of a uniform field.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::UniformField, point, quantity::VectorField; parameter::Parameter=Parameter())

a = excitation.amplitude

p = excitation.direction

return a * p
end

function field(excitation::UniformField, point, quantity::ScalarField; parameter::Parameter=Parameter())

    return -excitation.amplitude*dot(excitation.direction, point)
end
