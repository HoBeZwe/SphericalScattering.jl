
"""
    field(excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

Compute the field or potential of a uniform field.
"""
function field(excitation::UniformField, quantity::Field; parameter::Parameter=Parameter(), zeroRadius=0.0)

    F = zeros(fieldType(quantity), size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(excitation, point, quantity; parameter=parameter)
    end

    return F
end



"""
    field(excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field of a uniform field.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

    a = excitation.amplitude
    p = excitation.direction

    return a * p
end


"""
    field(excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential of a uniform field.
"""
function field(excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

    return -excitation.amplitude * dot(excitation.direction, point)
end
