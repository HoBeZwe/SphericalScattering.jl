
struct PlaneWave{T,C} <: Excitation 
    embedding::Medium{T}
    wavenumber::C
    amplitude::C
    direction::SVector{3,T}
    polarization::SVector{3,T}
end

planeWave(;
embedding    = Medium(ε0, μ0),
wavenumber   = error("missing argument `wavenumber`"),
amplitude    = 1.0,
direction    = SVector(0.0,0.0,-1.0),
polarization = SVector(1.0,0.0,0.0)
) = PlaneWave(embedding, wavenumber, amplitude, direction, polarization)