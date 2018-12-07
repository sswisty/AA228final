include("OrbitalMechFcns.jl")
include("MDPSAT_Helper.jl")


Rₑ = 6378.137
# Define the orbital elements of the satellite
a = Rₑ + 500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(0)
ω = 0
ν = deg2rad(87)
oe_sat = OrbitalElementVec(a,e,i,Ω,ω,ν)

# Define the orbital elements of the debris
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(90)
ω = 0
ν = deg2rad(87)
oe_deb = OrbitalElementVec(a,e,i,Ω,ω,ν)
