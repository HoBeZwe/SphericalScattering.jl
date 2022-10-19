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

function scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::VectorField; parameter::Parameter=Parameter())
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    E0 = excitation.amplitude

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding,sphere.filling),1)

    Bk = zeros(SMatrix{2,2,Float64},n)
    B = Matrix(I, 2,2)

    for k in range(n,stop=1, step = -1)
        B11=perms[k+1]+2*perms[k]
        B12=2*(perms[k+1]-perms[k])*a[k]^-3
        B21=(perms[k+1]-perms[k])*a[k]^3
        B22=2*perms[k+1]+perms[k]

        Bk[k]=1/(3*perms[k])*([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(Float64, n+1)
    Dk = zeros(Float64, n+1)

    Ck[1]=1
    Dk[n+1]=0
    Dk[1]=B[2,1]/B[1,1]
    Ck[n+1]=1/B[1,1]

    for k in range(n, stop=2, step=-1)
        Bk[k]
        [Ck[k+1],Dk[k+1]]
        Ck[k],Dk[k] = Bk[k] * [Ck[k+1],Dk[k+1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos+1]
    D = Dk[pos+1]

    return C * E0 *dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 .- E0
end

function scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::ScalarField; parameter::Parameter=Parameter())
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding,sphere.filling),1)

    Bk = zeros(SMatrix{2,2,Float64},n)
    B = Matrix(I, 2,2)

    for k in range(n,stop=1, step = -1)
        B11=perms[k+1]+2*perms[k]
        B12=2*(perms[k+1]-perms[k])*a[k]^-3
        B21=(perms[k+1]-perms[k])*a[k]^3
        B22=2*perms[k+1]+perms[k]

        Bk[k]=1/(3*perms[k])*([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(Float64, n+1)
    Dk = zeros(Float64, n+1)

    Ck[1]=1
    Dk[n+1]=0
    Dk[1]=B[2,1]/B[1,1]
    Ck[n+1]=1/B[1,1]

    for k in range(n, stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k],Dk[k] = Bk[k] * [Ck[k+1],Dk[k+1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos+1]
    D = Dk[pos+1]
    return C * Φ0 - D * Φ0  / r^3 - Φ0
end

function scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField, point, quantity::ScalarField; parameter::Parameter=Parameter())
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    perms = getfield.(vcat(sphere.embedding,sphere.filling),1)

    Bk = zeros(SMatrix{2,2,Float64},n)
    B = Matrix(I, 2,2)

    for k in range(n,stop=1, step = -1)
        B11=perms[k+1]+2*perms[k]
        B12=2*(perms[k+1]-perms[k])*a[k]^-3
        B21=(perms[k+1]-perms[k])*a[k]^3
        B22=2*perms[k+1]+perms[k]

        Bk[k]=1/(3*perms[k])*([B11 B12; B21 B22])
        B = Bk[k] * B
    end
    Ck = zeros(Float64, n+1)
    Dk = zeros(Float64, n+1)

    Ck[1]=1
    Dk[1]=(B[2,1]+B[2,2]*a[end]^3)/(B[1,1]+B[1,2]*a[end]^3)
    Ck[n+1]=1/(B[1,1]+B[1,2]*a[end]^3)
    Dk[n+1]=a[end]^3*Ck[n+1]

    for k in range(n, stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k],Dk[k] = Bk[k] * [Ck[k+1],Dk[k+1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:n+1
        if r < a[k]
            pos = k
        end
    end

    if pos == n+1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos+1]
        D = Dk[pos+1]
        R = a[pos+1]
    end

    return C*Φ0- D * Φ0/ r^3 - Φ0

end

function scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField, point, quantity::VectorField; parameter::Parameter=Parameter())
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    E0 = excitation.amplitude
    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)-1
    perms = getfield.(vcat(sphere.embedding,sphere.filling),1)

    Bk = zeros(SMatrix{2,2,Float64},n)
    B = Matrix(I, 2,2)

    for k in range(n,stop=1, step = -1)
        B11=perms[k+1]+2*perms[k]
        B12=2*(perms[k+1]-perms[k])*a[k]^-3
        B21=(perms[k+1]-perms[k])*a[k]^3
        B22=2*perms[k+1]+perms[k]

        Bk[k]=1/(3*perms[k])*([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(Float64, n+1)
    Dk = zeros(Float64, n+1)

    Ck[1]=1
    Dk[1]=(B[2,1]+B[2,2]*a[end]^3)/(B[1,1]+B[1,2]*a[end]^3)
    Ck[n+1]=1/(B[1,1]+B[1,2]*a[end]^3)
    Dk[n+1]=a[end]^3*Ck[n+1]

    for k in range(n, stop=2, step=-1)
        Ck[k],Dk[k] = Bk[k] * [Ck[k+1],Dk[k+1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:n+1
        if r < a[k]
            pos = k
        end
    end

    if pos == n+1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos+1]
        D = Dk[pos+1]
        R = a[pos+1]
    end

    return C *E0*dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 - E0*dir

end

fieldType(F::VectorField) = SVector{3,Complex{eltype(F.locations[1])}}
fieldType(F::ScalarField) = Complex{eltype(F.locations[1])}
