"""
Structure for easy OE implementation
    define all in Radians always!!!
"""

struct OrbitalElementVec
    # Semi-major axis
    a::Float64
    # Eccentricity
    e::Float64
    # Inclination
    i::Float64
    # Right Ascension of Ascending Node
    Ω::Float64
    # Argument of Periapsis
    ω::Float64
    # True Anomaly
    ν::Float64
end
