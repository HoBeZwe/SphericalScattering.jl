
"""
    SoundHard

The normal velocity, and hence the normal derivative of the total pressure, vanishes on the surface.
"""
struct SoundHard <: AcousticBoundary end

"""
    SoundSoft

The total pressure vanishes on the surface (pressure release).
"""
struct SoundSoft <: AcousticBoundary end
