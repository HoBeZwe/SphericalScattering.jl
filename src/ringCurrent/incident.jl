

"""
    field(excitation::ElectricRingCurrent, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field radiated by an electric ring current at some position and orientation
"""
function field(excitation::ElectricRingCurrent, quantity::Field; parameter::Parameter=Parameter())

    F = zeros(SVector{3,Complex{Float64}}, length(quantity.locations))

    points = translate(quantity.locations, excitation.center)
    rotate!(points, -excitation.rotation)

    for (ind, point) in enumerate(points)
        F[ind] = field(excitation, point, quantity, parameter=parameter)
    end

    rotate!(F, excitation.rotation)

    return F
end



"""
    field(excitation::ElectricRingCurrent, r, ϑ, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field radiated by an electric ring current placed in origin at point (r, ϑ) 
"""
function field(excitation::ElectricRingCurrent, point, quantity::ElectricField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = excitation.radius

    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε

    eps = parameter.relativeAccuracy
    
    Eϕ = complex(0.0) # initialize
    δE = Inf
    n  = -1
     
    r = point_sph[1]

    kr = k * r 
    kR = k * R

    kr == 0.0 && return SVector(complex(0.0), complex(0.0), complex(0.0)) # in origin / at f=0 Hz the electric field is zero
    
    cosθ = 0.0
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    

    try
        if r < R                
            while δE > eps || n < 10
                n += 2                      # (for z0=0 the even expansion_coeff are 0)
                expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                
                Jkr = besselj(n + 0.5, kr)  # Bessel function 1st kind 
                HkR = hankelh2(n + 0.5, kR) # Hankel function 2nd kind
                
                ΔE = expansion_coeff * HkR * Jkr * dnPl(cosϑ, n, 1)
                δE = abs(ΔE / Eϕ)           # relative change 

                Eϕ += ΔE 
            end
        
        else # r >=R
            while δE > eps || n < 10
                n += 2                      # (for z0=0 the even expansion_coeff are 0)
                expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                                        
                Hkr = hankelh2(n + 0.5, kr) # Hankel function 2nd kind
                JkR = besselj(n + 0.5, kR)  # Bessel function 1st kind     
                        
                ΔE = expansion_coeff * Hkr * JkR * dnPl(cosϑ, n, 1)
                δE = abs(ΔE / Eϕ)           # relative change 

                Eϕ += ΔE
            end
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end 

    Eϕ *= I0 * sqrt(μ / ε) * sinϑ * k * sqrt(R / r) / 2 * π # constant factors
    
    #return SVector(0.0,0.0,Eϕ)
    return SVector(-Eϕ*sin(point_sph[3]), Eϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end



"""
    field(excitation::ElectricRingCurrent, r, ϑ, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field radiated by an electric ring current placed in origin at point (r, ϑ) 
"""
function field(excitation::ElectricRingCurrent, point, quantity::MagneticField; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = excitation.radius

    eps = parameter.relativeAccuracy
    
    Hr  = complex(0.0) # initialize
    Hϑ  = complex(0.0) # initialize
    δHr = Inf
    n   = -1
     
    r = point_sph[1]

    kr = k * r 
    kR = k * R

    kr == 0.0 && return SVector(complex(0.0), complex(0.0), complex(0.0)) # in origin / at f=0 Hz the field is zero
    
    cosθ = 0.0
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    

    try
        if r <= R                
            while δHr > eps || n < 10
                n += 2                      # (for z0=0 the even expansion_coeff are 0)
                expansion_coeff = (2 * n + 1) / 2  * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                
                jkr = besselj(n + 0.5, kr)   # Bessel function 1st kind (without sqrt-term)
                hkR = hankelh2(n + 0.5, kR)  # Hankel function 2nd kind (without sqrt-term)
                jkrt = besselj(n + 1.5, kr)  # for derivative of Riccati-Bessel function 1st kind

                dJkr = (n + 1) * sqrt(R / r) * jkr * hkR - k * sqrt(r * R) * hkR * jkrt  # HkR * derivative of Riccati-Bessel function 1st kind
                
                ΔHr = expansion_coeff * hkR * jkr * Pl(cosϑ, n)
                ΔHϑ = expansion_coeff / n / (n + 1) * dJkr * dnPl(cosϑ, n, 1)

                δHr = abs(ΔHr / Hr)           # relative change 

                Hr += ΔHr
                Hϑ += ΔHϑ 
            end
        
        else # r >=R
            while δHr > eps || n < 10
                n += 2                       # (for z0=0 the even expansion_coeff are 0)
                expansion_coeff = (2 * n + 1) / 2 * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                                        
                hkr = hankelh2(n + 0.5, kr)  # Hankel function 2nd kind (without sqrt-term)
                jkR = besselj(n + 0.5, kR)   # Bessel function 1st kind (without sqrt-term)
                hkrt = hankelh2(n + 1.5, kr) # for derivative of Riccati-Hankel function 2nd kind
                
                dHkr = (n + 1) * sqrt(R / r) * hkr * jkR - k * sqrt(r * R) * jkR * hkrt  # JkR * derivative of Riccati-Hankel function 2nd kind
                        
                ΔHr = expansion_coeff * hkr * jkR * Pl(cosϑ, n)
                ΔHϑ = expansion_coeff / n / (n + 1) * dHkr * dnPl(cosϑ, n, 1)
     
                δHr = abs(ΔHr / Hr)           # relative change

                Hr += ΔHr
                Hϑ += ΔHϑ
            end
        end
    catch
        # print("did not converge: n=$n\n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (not verified yet)
    end
    
    Hr *= im * I0 * sqrt(R / r) / r * π / 2
    Hϑ *= -im * I0 * sinϑ / r * π / 2

    return convertSpherical2Cartesian(SVector(Hr, Hϑ, 0.0), point_sph)
end



"""
    field(excitation::ElectricRingCurrent, r, ϑ, quantity::FarField; parameter::Parameter=Parameter())

Compute the (electric) far-field radiated by an electric ring current placed in origin at point (r, ϑ) 
"""
function field(excitation::ElectricRingCurrent, point, quantity::FarField; parameter::Parameter=Parameter())
    
    point_sph = cart2sph(point) # [r ϑ φ]

    k  = excitation.wavenumber
    I0 = excitation.amplitude
    R  = excitation.radius
    
    μ  = excitation.embedding.μ
    ε  = excitation.embedding.ε
    
    eps = parameter.relativeAccuracy
    
    Eϕ = complex(0.0) # initialize
    δE = Inf
    n  = -1
    
    r = point_sph[1]

    kr = k * r 
    kR = k * R
    
    cosθ = 0.0
    sinϑ = sin(point_sph[2]) 
    cosϑ = cos(point_sph[2])    
    
    try
        while δE > eps || n < 10
            n += 2                      # (for z0=0 the even expansion_coeff are 0)
            expansion_coeff = (2 * n + 1) / 2 / n / (n + 1) * dnPl(cosθ, n, 1) # variable part of expansion coefficients
                                        
            Hkr = hankelh2(n + 0.5, kr) # Hankel function 2nd kind
            JkR = besselj(n + 0.5, kR)  # Bessel function 1st kind     
                        
            ΔE = expansion_coeff * Hkr * JkR * dnPl(cosϑ, n, 1)
            δE = abs(ΔE / Eϕ)           # relative change 

            Eϕ += ΔE
        end
    catch
        # print("did not converge: n=$n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified)
    end 
    
    Eϕ *= -I0 * sqrt(μ / ε) * sinϑ * k * sqrt(R / r) / 2 * π # constant factors
    
    return SVector(-Eϕ*sin(point_sph[3]), Eϕ*cos(point_sph[3]), 0.0) # convert to Cartesian representation
end