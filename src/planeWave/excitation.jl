
struct PlaneWave{T,R,C} <: Excitation
    embedding::Medium{C}
    wavenumber::R
    amplitude::T
    direction::SVector{3,R}
    polarization::SVector{3,R}
end

planeWave(;
    embedding    = Medium(ε0, μ0),
    wavenumber   = error("missing argument `wavenumber`"),
    amplitude    = 1.0,
    direction    = SVector{3,typeof(wavenumber)}(0.0, 0.0, -1.0),
    polarization = SVector{3,typeof(wavenumber)}(1.0, 0.0, 0.0),
) = PlaneWave(embedding, wavenumber, amplitude, direction, polarization)
