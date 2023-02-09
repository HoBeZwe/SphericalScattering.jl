
"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a sphere, for an incident uniform field.
"""
function scatteredfield(sphere::Sphere, excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    F = zeros(fieldType(quantity), size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, point, quantity; parameter=parameter)
    end

    return F
end



"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a Dielectric sphere, for an incident uniform field with polarization in given direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter()
)

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return (3 * ε0 / (ε1 + 2 * ε0) - 1) * E0 # Sihvola&Lindell 1988, (8)
    end

    return E0 * (-(ε1 - ε0) / (ε1 + 2 * ε0) * R^3 / r^3) + 3 * (ε1 - ε0) / (ε1 + 2 * ε0) * R^3 * point / r^5 * dot(E0, point) #Sihvola&Lindell 1988, (9)
end



"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a Dielectric sphere, for an incident uniform field with polarization in given direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter()
)

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return (3 * ε0 / (ε1 + 2 * ε0) - 1) * Φ0 # Sihvola&Lindell 1988, (8)
    end

    return -Φ0 * ((ε1 - ε0) / (ε1 + 2 * ε0) * R^3 / r^3) #Sihvola&Lindell 1988, (9)
end



"""
    scatteredfield(sphere::DielectricSphereThinImpedanceLayer, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a dielectric sphere with a thin coating, where the displacement field in the coating is only in radial direction.
We assume an an incident uniform field with polarization in the given direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinImpedanceLayer,
    excitation::UniformField,
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
)

    point = rotate(excitation, [point]; inverse=true)[1]
    point_sph = cart2sph(point)

    r = norm(point)

    ẑ = SVector(0.0, 0.0, 1.0) #excitation.direction
    cosθ = dot(ẑ, point) / r
    sinθ = norm(cross(ẑ, point)) / r

    A, K = scatterCoeff(sphere, excitation)

    if r > sphere.radius
        E = -SVector((-2 * A / r^3) * cosθ, (+A / r^3) * (-sinθ), 0.0)
    else
        E = -SVector((-K * cosθ, K * sinθ, 0.0))
    end

    E_cart = convertSpherical2Cartesian(E, point_sph)

    return rotate(excitation, [E_cart]; inverse=false)[1]
end



"""
    scatteredfield(sphere::DielectricSphereThinImpedanceLayer, excitation::UniformField, point, quantity::DisplacementField; parameter::Parameter=Parameter())

Compute the displacement field D = ε * E.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinImpedanceLayer,
    excitation::UniformField,
    point,
    quantity::DisplacementField;
    parameter::Parameter=Parameter(),
)

    E = scatteredfield(sphere, excitation, point, ElectricField(quantity.locations); parameter=parameter)

    if norm(point) > sphere.radius
        D = sphere.embedding.ε * E
    else
        D = sphere.filling.ε * E
    end

    return D
end



"""
    scatteredfield(sphere::DielectricSphereThinImpedanceLayer, excitation::UniformField, point, quantity::ScalarPotentialJump; parameter::Parameter=Parameter())

Compute the jump of the scalar potential for a dielectric sphere with a thin coating, where the displacement field in the coating is only in radial direction.
We assume an an incident uniform field with polarization in the given direction.

More precisely, we compute the difference Δ = Φ_i - Φ_e, where Φ_i is the potential on the inner side and ϕ_e the exterior potential.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinImpedanceLayer,
    excitation::UniformField,
    point,
    quantity::ScalarPotentialJump;
    parameter::Parameter=Parameter(),
)

    cosθ = dot(excitation.direction, point) / norm(point)

    ~, K = scatterCoeff(sphere, excitation)

    return sphere.thickness * (sphere.filling.ε / sphere.thinlayer.ε) * K * cosθ

end



"""
    scatteredfield(sphere::DielectricSphereThinImpedanceLayer, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a dielectric sphere with a thin coating, where the displacement field in the coating is only in radial direction.
We assume an an incident uniform field with polarization in the given direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinImpedanceLayer,
    excitation::UniformField,
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
)

    cosθ = dot(excitation.direction, point) / norm(point)
    r = norm(point)

    R = sphere.radius

    A, K = scatterCoeff(sphere, excitation)

    if r >= R
        return (+A / r^2) * cosθ # Jones 1995, (C.1a)
    else
        return (-K * r * cosθ)
    end
end



"""
    scatterCoeff(sp::DielectricSphereThinImpedanceLayer, ex::UniformField)

Compute the expansion coefficients for the thin impedance layer case and a uniform static field excitation.
"""
function scatterCoeff(sp::DielectricSphereThinImpedanceLayer, ex::UniformField)
    R = sp.radius
    Δ = sp.thickness
    εₘ = sp.thinlayer.ε
    εₑ = sp.embedding.ε
    εᵢ = sp.filling.ε
    E₀ = ex.amplitude

    A = E₀ * R^3 * (-R * εₑ * εₘ + R * εᵢ * εₘ - Δ * εₑ * εᵢ) / (2R * εₑ * εₘ + R * εᵢ * εₘ + 2Δ * εₑ * εᵢ)
    K = 3E₀ * R * εₑ * εₘ / (2 * R * εₘ * εₑ + R * εᵢ * εₘ + 2Δ * εₑ * εᵢ)

    return A, K
end



"""
    scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, for an incident uniform field with polarization in the given direction.

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return -E0
    end

    return -E0 * R^3 / r^3 + 3 * R^3 / r^5 * point * dot(E0, point) # Griffits, Example 3.8
end



"""
    scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a PEC sphere, for an incident uniform field with polarization in the given direction.

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::PECSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter()
)

    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return -Φ0
    end

    return Φ0 * (-R^3 / r^3) # Griffits, Example 3.8
end



"""
    scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a layered dielectric sphere, for an incident uniform field with polarization in the given direction
using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`.

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSphere{LN,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
) where {LN,LR,LC,FC,FT,FR}

    E0 = excitation.amplitude

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[n + 1] = 0
    Dk[1] = B[2, 1] / B[1, 1]
    Ck[n + 1] = 1 / B[1, 1]

    for k in range(n; stop=2, step=-1)
        Bk[k]
        [Ck[k + 1], Dk[k + 1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos + 1]
    D = Dk[pos + 1]

    return C * E0 * dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 - E0 * dir
end



"""
    scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a layered dielectric sphere, for an incident uniform field with polarization in the given direction
using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`.

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSphere{LN,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
) where {LN,LR,LC,FC,FT,FR}

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[n + 1] = 0
    Dk[1] = B[2, 1] / B[1, 1]
    Ck[n + 1] = 1 / B[1, 1]

    for k in range(n; stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos + 1]
    D = Dk[pos + 1]

    return C * Φ0 - D * Φ0 / r^3 - Φ0
end



"""
    scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField{FC,FT,FR}, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a layered dielectric sphere with PEC core, for an incident uniform field with polarization in the given direction
using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSpherePEC{LN,LD,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
) where {LN,LD,LR,LC,FC,FT,FR}

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a) - 1
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end
    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[1] = (B[2, 1] + B[2, 2] * a[end]^3) / (B[1, 1] + B[1, 2] * a[end]^3)
    Ck[n + 1] = 1 / (B[1, 1] + B[1, 2] * a[end]^3)
    Dk[n + 1] = a[end]^3 * Ck[n + 1]

    for k in range(n; stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:(n + 1)
        if r < a[k]
            pos = k
        end
    end

    if pos == n + 1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos + 1]
        D = Dk[pos + 1]
        R = a[pos + 1]
    end

    return C * Φ0 - D * Φ0 / r^3 - Φ0

end



"""
    scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a layered dielectric sphere with PEC core, for an incident uniform field with polarization in the given direction.

The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSpherePEC{LN,LD,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
) where {LN,LD,LR,LC,FC,FT,FR}
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    E0 = excitation.amplitude
    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a) - 1
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[1] = (B[2, 1] + B[2, 2] * a[end]^3) / (B[1, 1] + B[1, 2] * a[end]^3)
    Ck[n + 1] = 1 / (B[1, 1] + B[1, 2] * a[end]^3)
    Dk[n + 1] = a[end]^3 * Ck[n + 1]

    for k in range(n; stop=2, step=-1)
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:(n + 1)
        if r < a[k]
            pos = k
        end
    end

    if pos == n + 1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos + 1]
        D = Dk[pos + 1]
        R = a[pos + 1]
    end

    return C * E0 * dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 - E0 * dir

end

fieldType(F::ElectricField)       = SVector{3,Complex{eltype(F.locations[1])}}
fieldType(F::DisplacementField)   = SVector{3,Complex{eltype(F.locations[1])}}
fieldType(F::ScalarPotential)     = Complex{eltype(F.locations[1])}
fieldType(F::ScalarPotentialJump) = Complex{eltype(F.locations[1])}
