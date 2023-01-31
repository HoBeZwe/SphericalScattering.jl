
"""
    field(excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field radiated by a magnetic/electric ring current at some position and with some orientation.
"""
function field(excitation::Dipole, quantity::Field; parameter::Parameter=Parameter())

    T = typeof(excitation.frequency)

    F = zeros(SVector{3,Complex{T}}, length(quantity.locations))

    # --- distinguish electric/magnetic current
    fieldType, exc = getFieldType(excitation, quantity)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = field(exc, point, fieldType; parameter=parameter)
    end

    return F
end



"""
    field(excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field radiated by a Hertzian dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::Dipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = wavenumber(excitation)
    ε  = excitation.embedding.ε
    μ  = excitation.embedding.μ

    r0 = excitation.position
    p  = excitation.orientation

    d = point - r0
    r = norm(d)
    n = d / r

    # r < 0.9 && return SVector(NaN, NaN, NaN)

    return Il / (4 * π) * sqrt(μ / ε) * exp(-im * k * r) * (k / r * ((n × p) × n) + (1 / k / r^3 + im / r^2) * (3 * n * dot(n, p) - p))
end



"""
    field(excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

Compute the magnetic field radiated by a Hertzian dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::Dipole, point, quantity::MagneticField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = wavenumber(excitation)

    r0 = excitation.position
    p  = excitation.orientation

    d = point - r0
    r = norm(d)
    n = d / r

    return Il / (4 * π) * (n × p) * exp(-im * k * r) / r * (k + 1 / (im * r))
end



"""
    field(excitation::HertzianDipole, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far field radiated by a Hertzian dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::HertzianDipole, point, quantity::FarField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = wavenumber(excitation)
    ε  = excitation.embedding.ε
    μ  = excitation.embedding.μ

    r0 = excitation.position
    p  = excitation.orientation

    n = point / norm(point)
    kd = k * r0 ⋅ point         # phase correction if dipole is not placed in origin

    return Il / (4 * π) * sqrt(μ / ε) * k * ((n × p) × n) * exp(im * kd)
end



"""
    field(excitation::FitzgeraldDipole, point, quantity::FarField; parameter::Parameter=Parameter())

Compute the electric far field radiated by a Fitzgerald dipole at given position and orientation at point.

The point and the returned field are in Cartesian coordinates.
"""
function field(excitation::FitzgeraldDipole, point, quantity::FarField; parameter::Parameter=Parameter())

    Il = excitation.amplitude
    k  = wavenumber(excitation)

    r0 = excitation.position
    p  = excitation.orientation

    n = point / norm(point)
    kd = k * r0 ⋅ point         # phase correction if dipole is not placed in origin

    return Il / (4 * π) * (n × p) * k * exp(im * kd)
end



# function field(excitation::HertzianDipole, point, quantity::ElectricField; parameter::Parameter=Parameter())

#     point_sph = cart2sph(point) # [r ϑ φ]

#     k  = wavenumber(excitation)
#     I0 = excitation.amplitude
#     R  = norm(excitation.position)

#     μ  = excitation.embedding.μ
#     ε  = excitation.embedding.ε

#     eps = parameter.relativeAccuracy

#     Hr  = complex(0.0) # initialize
#     Hϑ  = complex(0.0) # initialize
#     δHr = Inf
#     n   = 0

#     r = point_sph[1]

#     kr = k * r 
#     kR = k * R

#     kr == 0.0 && return SVector(complex(0.0), complex(0.0), complex(0.0)) # in origin / at f=0 Hz the field is zero

#     sinϑ = sin(point_sph[2]) 
#     cosϑ = cos(point_sph[2])    

#     try
#         if r <= R                
#             while δHr > eps || n < 10
#                 n += 1
#                 expansion_coeff = 2 * n + 1  # variable part of expansion coefficients

#                 jkr = besselj(n + 0.5, kr)   # Bessel function 1st kind (without sqrt-term)
#                 hkR = hankelh2(n + 0.5, kR)  # Hankel function 2nd kind (without sqrt-term)
#                 jkrt = besselj(n + 1.5, kr)  # for derivative of Riccati-Bessel function 1st kind

#                 dJkr = (n + 1) * sqrt(R / r) * jkr * hkR - k * sqrt(r * R) * hkR * jkrt  # HkR * derivative of Riccati-Bessel function 1st kind

#                 ΔHr = expansion_coeff * n * (n + 1) * hkR * jkr * Pl(cosϑ, n)
#                 ΔHϑ = expansion_coeff * dJkr * dnPl(cosϑ, n, 1)

#                 δHr = abs(ΔHr / Hr)           # relative change 

#                 Hr += ΔHr
#                 Hϑ += ΔHϑ 
#                 # @show Hr
#             end

#         else # r >=R
#             while δHr > eps || n < 10
#                 n += 1
#                 expansion_coeff = 2 * n + 1  # variable part of expansion coefficients

#                 hkr = hankelh2(n + 0.5, kr)  # Hankel function 2nd kind (without sqrt-term)
#                 jkR = besselj(n + 0.5, kR)   # Bessel function 1st kind (without sqrt-term)
#                 hkrt = hankelh2(n + 1.5, kr) # for derivative of Riccati-Hankel function 2nd kind

#                 dHkr = (n + 1) * sqrt(R / r) * hkr * jkR - k * sqrt(r * R) * jkR * hkrt  # JkR * derivative of Riccati-Hankel function 2nd kind

#                 ΔHr = expansion_coeff * n * (n + 1) * hkr * jkR * Pl(cosϑ, n)
#                 ΔHϑ = expansion_coeff * dHkr * dnPl(cosϑ, n, 1)

#                 δHr = abs(ΔHr / Hr)           # relative change

#                 Hr += ΔHr
#                 Hϑ += ΔHϑ
#                 # @show Hr
#             end
#         end
#     catch
#         # print("did not converge: n=$n\n") # if Hankel function throws overflow error -> result still agrees with small argument approximation (verified yet)
#     end
#     # @show n
#     ZF = sqrt(μ / ε)
#     Hr *= -im*I0 * sqrt(R / r) / r / 8 / k / R / R * ZF
#     Hϑ *=  im*I0 * sinϑ / r / 8 / k / R / R * ZF

#     return convertSpherical2Cartesian(SVector(Hr, Hϑ, 0.0), point_sph)
# end
