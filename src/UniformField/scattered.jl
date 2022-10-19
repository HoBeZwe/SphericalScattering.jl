"""
    scatteredfield(sphere::PECSphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a sphere, for an incident uniform field.
"""
function scatteredfield(sphere::Sphere, excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    F = zeros(fieldType(quantity), size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, point, quantity, parameter=parameter)
    end

    return F
end



"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the vector field scattered by a Dielectric sphere, for an incident uniform field with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::VectorField; parameter::Parameter=Parameter())

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return  (3 * ε0 / (ε1 + 2 * ε0) -1) * E0 # Sihvola&Lindell 1988, (8)
    end

    return E0 * (-(ε1 -ε0)/(ε1 + 2*ε0) * R^3 / r^3) + 3 * (ε1 -ε0)/(ε1 + 2*ε0) * R^3 * point / r^5 * dot(E0,point) #Sihvola&Lindell 1988, (9)

end

"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the scalar field scattered by a Dielectric sphere, for an incident uniform field with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ScalarField; parameter::Parameter=Parameter())

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return  (3 * ε0 / (ε1 + 2 * ε0) -1)* Φ0 # Sihvola&Lindell 1988, (8)
    end

    return -Φ0 * ((ε1 -ε0)/(ε1 + 2*ε0) * R^3 / r^3) #Sihvola&Lindell 1988, (9)

end

function scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::VectorField; parameter::Parameter=Parameter())

    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return SVector{3, T}(0.0,0.0,0.0)-E0
    end

    return -E0*R^3/r^3 + 3*R^3/r^5*point*dot(E0,point) # Griffits, Example 3.8
end

function scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ScalarField; parameter::Parameter=Parameter())

    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return -Φ0
    end

    return Φ0 *(- R^3/r^3) # Griffits, Example 3.8
end

