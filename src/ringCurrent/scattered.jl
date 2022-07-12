

"""
    field(excitation::ElectricRingCurrent, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field radiated by an electric ring current at some position and orientation
"""
function scatteredfield(sphere::PECSphere, excitation::RingCurrent, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    T = typeof(excitation.wavenumber)

    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- translate/rotate coordinates
    points = quantity.locations #translate(quantity.locations, -excitation.center)
    #rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(points)
        F[ind] = scatteredfield(sphere, exc, point, fieldType, parameter=parameter)
    end

    # --- rotate resulting field
    #rotate!(F, excitation.rotation)

    return F
end



"""
    scatteredfield(sphere::PECSphere, excitation::RingCurrent, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by the PEC sphere, where the ring current is placed along the z-axis.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::RingCurrent, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = sqrt(norm(excitation.center)^2 + excitation.radius^2)    # distance loop-origin
    θ  = atan(excitation.radius / norm(excitation.center))        # angle between z-axis and connection loop-origin

    T = typeof(k)

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Eϕ = Complex{T}(0.0) # initialize
    δE  = T(Inf)
    n   = 0

    r = point_sph[1]
     
    kr = k * r
    kR = k * R
    ka = k * sphere.radius

    r <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0
    
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 1                      # (for z0=0 the even expansion_coeff are 0)
            expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
            
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            HkR = hankelh2(n + T(0.5), kR) # Hankel function 2nd kind
            Hkr = hankelh2(n + T(0.5), kr) # Hankel function 2nd kind
                    
            ΔE = expansion_coeff * Jka_Hka * HkR * Hkr * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δE = abs(ΔE / Eϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)
            
            Eϕ += ΔE 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϕ *= I0 * sqrt(μ / ε) * k / sqrt(R * r) * sinϑ * sinθ / 2 * π * excitation.radius

    return SVector{3,Complex{T}}(-Eϕ*sin(point_sph[3]), Eϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::RingCurrent, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field scattered by the PEC sphere, where the ring current is placed along the z-axis.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::RingCurrent, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = sqrt(norm(excitation.center)^2 + excitation.radius^2)    # distance loop-origin
    θ  = atan(excitation.radius / norm(excitation.center))        # angle between z-axis and connection loop-origin

    T = typeof(k)

    eps = parameter.relativeAccuracy
    
    Hr  = Complex{T}(0.0) # initialize
    Hϑ  = Complex{T}(0.0) # initialize
    δHr = T(Inf)
    n   = 0
    
    r = point_sph[1]

    kr = k * r 
    kR = k * R
    ka = k * sphere.radius

    r <= sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0
    
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    

    try
        while δHr > eps || n < 10
            n += 1                       # (for z0=0 the even expansion_coeff are 0)
            expansion_coeff = (2 * n + 1) / 2 * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                           
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            Hkr = hankelh2(n + T(0.5), kr)  # Hankel function 2nd kind (without sqrt-term)
            HkR = hankelh2(n + T(0.5), kR)  # Hankel function 2nd kind (without sqrt-term)
            Hkrt = hankelh2(n + T(1.5), kr) # for derivative of Riccati-Hankel function 2nd kind
  
            dHkr = (n + 1) * sqrt(R / r) * Jka_Hka * Hkr * HkR - k * sqrt(r * R) * Jka_Hka * HkR * Hkrt  # HkR * derivative of Riccati-Hankel function 2nd kind
                        
            ΔHr = expansion_coeff * Jka_Hka * Hkr * HkR * Pl(cosϑ, n)
            ΔHϑ = expansion_coeff / n / (n + 1) * dHkr * dnPl(cosϑ, n, 1)
     
            isodd(n) &&  (δHr = abs(ΔHr / Hr)) # relative change every second n (for z0=0 the even expansion_coeff are 0)
            
            Hr += ΔHr
            Hϑ += ΔHϑ
        end
    catch
        # print("did not converge: n=$n\n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Hr *= -im * I0 * excitation.radius * sinθ * sqrt(r / R) / r / r * π / 2 #* μ / ε
    Hϑ *=  im * I0 * excitation.radius * sinϑ * sinθ / r / R * π / 2 #* μ / ε

    return convertSpherical2Cartesian(SVector{3,Complex{T}}(Hr, Hϑ, 0.0), point_sph)
end



"""
    scatteredfield(sphere::PECSphere, excitation::ElectricRingCurrent, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by the PEC sphere, where the ring current is placed along the z-axis.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::ElectricRingCurrent, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]
    
    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = sqrt(norm(excitation.center)^2 + excitation.radius^2)    # distance loop-origin
    θ  = atan(excitation.radius / norm(excitation.center))        # angle between z-axis and connection loop-origin

    T = typeof(k)

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Eϕ = Complex{T}(0.0) # initialize
    δE = T(Inf)
    n  = 0
     
    kR = k * R
    ka = k * sphere.radius
    
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 1                      # (for z0=0 the even expansion_coeff are 0)
            expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                    
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            HkR = hankelh2(n + T(0.5), kR) # Hankel function 2nd kind
                    
            ΔE = expansion_coeff * Jka_Hka * HkR * im^(n+1) * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δE = abs(ΔE / Eϕ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϕ += ΔE 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϕ *= I0 * sqrt(μ / ε) * excitation.radius * sinθ * sinϑ * k * sqrt(π / 2 / k / R)

    return SVector{3,Complex{T}}(-Eϕ*sin(point_sph[3]), Eϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    scatteredfield(sphere::PECSphere, excitation::MagneticRingCurrent, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far-field scattered by the PEC sphere, where the ring current is placed along the z-axis.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::MagneticRingCurrent, point, quantity::FarField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]
    
    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = sqrt(norm(excitation.center)^2 + excitation.radius^2)    # distance loop-origin
    θ  = atan(excitation.radius / norm(excitation.center))        # angle between z-axis and connection loop-origin

    T = typeof(k)

    eps = parameter.relativeAccuracy
    
    Eϑ = Complex{T}(0.0) # initialize
    δE = T(Inf)
    n  = 0
     
    kR = k * R
    ka = k * sphere.radius
    
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 1                      # (for z0=0 the even expansion_coeff are 0)
            expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                    
            Jka_Hka = scatterCoeff(sphere, excitation, n, ka)

            HkR = hankelh2(n + T(0.5), kR) # Hankel function 2nd kind
                    
            ΔE = expansion_coeff * Jka_Hka * HkR * im^(n+1) * dnPl(cosϑ, n, 1)
            isodd(n) &&  (δE = abs(ΔE / Eϑ)) # relative change every second n (for z0=0 the even expansion_coeff are 0)

            Eϑ += ΔE 
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end

    Eϑ *= I0 * excitation.radius * sinθ * sinϑ * k * sqrt(π / 2 / k / R)

    return SVector{3,Complex{T}}(Eϑ*cos(point_sph[2])*cos(point_sph[3]), Eϑ*cos(point_sph[2])*sin(point_sph[3]), -Eϑ*sin(point_sph[2])) # convert to Cartesian representation
end



"""
    scatterCoeff(sphere::PECSphere, excitation::ElectricRingCurrent, n::Int, ka)

Compute scattering coefficient for electric ring current.
"""
function scatterCoeff(sphere::PECSphere, excitation::ElectricRingCurrent, n::Int, ka)

    T = typeof(excitation.wavenumber)

    Jka = besselj(n + T(0.5), ka)  # Bessel function 1st kind
    Hka = hankelh2(n + T(0.5), ka) # Hankel function 2nd kind

    return Jka / Hka
end



"""
    scatterCoeff(sphere::PECSphere, excitation::MagneticRingCurrent, n::Int, ka)

Compute scattering coefficient for magnetic ring current.
"""
function scatterCoeff(sphere::PECSphere, excitation::MagneticRingCurrent, n::Int, ka)

    T = typeof(excitation.wavenumber)

    Jka = (n + 1) * besselj(n + T(0.5), ka)   - ka * besselj(n + T(1.5), ka)  # derivative spherical Bessel function 1st kind (without sqrt factor)
    Hka = (n + 1) * hankelh2(n + T(0.5), ka)  - ka * hankelh2(n + T(1.5), ka) # derivative spherical Hankel function 2nd kind (without sqrt factor)

    return Jka / Hka
end