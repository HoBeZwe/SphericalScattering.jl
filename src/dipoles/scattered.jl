
"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a dipole at some position and orientation.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    F = zeros(SVector{3,Complex{Float64}}, size(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- translate/rotate coordinates
    points = quantity.locations # translate(quantity.locations, -excitation.center)
    # rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(points)
        F[ind] = scatteredfield(sphere, exc, point, fieldType, parameter=parameter)
    end

    # --- rotate resulting field
    # rotate!(F, excitation.rotation)

    return F
end



"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    Il = excitation.amplitude
    z0 = norm(excitation.center)   # distance Dipole-origin

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Er  = complex(0.0) # initialize
    Eϑ  = complex(0.0) # initialize
    δEr = Inf
    n   = 0
    
    r = point_sph[1]

    kr  = k * r 
    kz0 = k * z0
    ka  = k * sphere.radius

    r <= sphere.radius && return SVector(complex(0.0), complex(0.0), complex(0.0)) # inside the sphere the field is 0
    
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    

    try
        while δEr > eps || n < 10
            n += 1                       
            expansion_coeff = 2 * n + 1  # variable part of expansion coefficients
                           
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkr = hankelh2(n + 0.5, kr)   # Hankel function 2nd kind (without sqrt-term)
            Hkz0 = hankelh2(n + 0.5, kz0) # Hankel function 2nd kind (without sqrt-term)
            Hkrt = hankelh2(n + 1.5, kr)  # for derivative of Riccati-Hankel function 2nd kind
  
            dHkr = (n + 1) * sqrt(z0 / r) * Jka_Hka * Hkr * Hkz0 - k * sqrt(r * z0) * Jka_Hka * Hkz0 * Hkrt  # HkR * derivative of Riccati-Hankel function 2nd kind
                        
            ΔEr = expansion_coeff * n * (n + 1) * Jka_Hka * Hkr * Hkz0 * Pl(cosϑ, n)
            ΔEϑ = expansion_coeff * dHkr * dnPl(cosϑ, n, 1)
     
            isodd(n) && (δEr = abs(ΔEr / Er)) # relative change every second n (for z0=0 the even expansion_coeff are 0)
            
            Er += ΔEr
            Eϑ += ΔEϑ
        end
    catch
        # print("did not converge: n=$n\n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    ZF = sqrt(μ / ε)

    Er *= -im * Il / 8 / k / z0 * sqrt(r / z0) / r / r * ZF  
    Eϑ *=  im * Il / 8 / k / z0 / z0 * sinϑ / r * ZF 

    return convertSpherical2Cartesian(SVector(Er, Eϑ, 0.0), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    Il = excitation.amplitude
    z0 = norm(excitation.center)    # distance Dipole-origin

    eps = parameter.relativeAccuracy
    
    Hϕ = complex(0.0) # initialize
    δH  = Inf
    n   = 0

    r = point_sph[1]
     
    kr  = k * r
    kz0 = k * z0
    ka  = k * sphere.radius

    r <= sphere.radius && return SVector(complex(0.0), complex(0.0), complex(0.0)) # inside the sphere the field is 0
    
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δH > eps || n < 10
            n += 1                      
            expansion_coeff = 2 * n + 1 # variable part of expansion coefficients
            
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + 0.5, kz0) # Hankel function 2nd kind
            Hkr  = hankelh2(n + 0.5, kr) # Hankel function 2nd kind
                    
            ΔH = expansion_coeff * Jka_Hka * Hkz0 * Hkr * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δH = abs(ΔH / Hϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)
            
            Hϕ += ΔH 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Hϕ *= -Il / z0 / sqrt(z0 * r) * sinϑ / 8

    return SVector(-Hϕ*sin(point_sph[3]), Hϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::HertzianDipole, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::HertzianDipole, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]
    
    k  = excitation.wavenumber
    Il = excitation.amplitude
    z0 = norm(excitation.center)    # distance dipole-origin

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Eϑ = complex(0.0) # initialize
    δE = Inf
    n  = 0
     
    kz0 = k * z0
    ka  = k * sphere.radius
    
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 1                      
            expansion_coeff = 2 * n + 1   # variable part of expansion coefficients
                    
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + 0.5, kz0) # Hankel function 2nd kind
                    
            ΔE = expansion_coeff * Jka_Hka * Hkz0 * im^(n+1) * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δE = abs(ΔE / Eϑ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϑ += ΔE 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϑ *= Il / k * sinϑ / z0 / 4 * sqrt(k / z0 / 2 / π) * sqrt(μ / ε)

    return SVector(Eϑ*cos(point_sph[2])*cos(point_sph[3]), Eϑ*cos(point_sph[2])*sin(point_sph[3]), -Eϑ*sin(point_sph[2])) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::FitzgeraldDipole, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by a PEC sphere, where the dipole is placed along the z-axis at z0.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::FitzgeraldDipole, point, quantity::FarField; parameter::Parameter=Parameter())
    
    point_sph = cart2sph(point) # [r ϑ φ]
    
    k  = excitation.wavenumber
    Ul = excitation.amplitude
    z0 = norm(excitation.center)    # distance dipole-origin

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Eϕ = complex(0.0) # initialize
    δE = Inf
    n  = 0
     
    kz0 = k * z0
    ka  = k * sphere.radius
    
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 1
            expansion_coeff = 2 * n + 1 # variable part of expansion coefficients
                    
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkz0 = hankelh2(n + 0.5, kz0) # Hankel function 2nd kind
                    
            ΔE = expansion_coeff * Jka_Hka * Hkz0 * im^(n+1) * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δE = abs(ΔE / Eϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϕ += ΔE 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϕ *= -Ul / k * sinϑ / z0 / 4 * sqrt(k / z0 / 2 / π)

    return SVector(-Eϕ*sin(point_sph[3]), Eϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatterCoeff(sphere::PECSphere, excitation::ElectricRingCurrent, n::Int, ka)

Compute scattering coefficient for electric ring current.
"""
function scatterCoeff(sphere::PECSphere, excitation::FitzgeraldDipole, n::Int, ka)

    Jka = besselj(n + 0.5, ka)  # Bessel function 1st kind 
    Hka = hankelh2(n + 0.5, ka) # Hankel function 2nd kind

    return Jka / Hka
end



"""
    scatterCoeff(sphere::PECSphere, excitation::MagneticRingCurrent, n::Int, ka)

Compute scattering coefficient for magnetic ring current.
"""
function scatterCoeff(sphere::PECSphere, excitation::HertzianDipole, n::Int, ka)

    Jka = (n + 1) * besselj(n + 0.5, ka)   - ka * besselj(n + 1.5, ka)  # derivative spherical Bessel function 1st kind (without sqrt factor)
    Hka = (n + 1) * hankelh2(n + 0.5, ka)  - ka * hankelh2(n + 1.5, ka) # derivative spherical Hankel function 2nd kind (without sqrt factor)

    return Jka / Hka
end