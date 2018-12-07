include("OrbitMechFcns.jl")
include("MDPSAT_Helper.jl")

Rₑ = 6378.137
μ = 3.986e5
# Define the orbital elements of the satellite
a = Rₑ + 500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(0)
Ω = deg2rad(0)
ω = 0
ν = deg2rad(359.5)
oe_sat = OrbitalElementVec(a,e,i,Ω,ω,ν)

# Define the orbital elements of the debris
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(0.1)
Ω = deg2rad(0)
ω = 0
ν = deg2rad(359.5)
oe_deb = OrbitalElementVec(a,e,i,Ω,ω,ν)

A = [-10.0, 0.0, 10.0]
s0 = initState(oe_sat, oe_deb)
mdp = SATMDP(-10, -1, Rₑ, μ, 0, 0.9, A, 1)
depth = 10

action_taken = []
relDist = []
si = []
rewards = []
push!(relDist, norm(s0[1:3] - s0[7:9]))
push!(si, s0)

for i = 1:7
    a, v = selectAction(si[i], depth,mdp)
    sp = simulateDynamics(si[i], a, mdp)
    rew = rewardModel(si[i], a, mdp)
    push!(action_taken, a)
    push!(relDist, norm(sp[1:3] - sp[7:9]))
    push!(rewards, rew)
    push!(si, sp)
end
