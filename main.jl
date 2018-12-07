include("OrbitMechFcns.jl")
include("MDPSAT_Helper.jl")

## Constants
Rₑ = 6378.137
μ = 3.986e5

## Orbit Initialization
h = 500 # altitude
# Orbital elements of the satellite
a = Rₑ + h; e = 0.0; i = deg2rad(0); Ω = deg2rad(0); ω = 0; ν = deg2rad(359.5)
oe_sat = OrbitalElementVec(a, e, i, Ω, ω, ν)
# Orbital elements of the debris
a = Rₑ + h; e = 0.0; i = deg2rad(0.1); Ω = deg2rad(0); ω = 0; ν = deg2rad(359.5)
oe_deb = OrbitalElementVec(a, e, i, Ω, ω, ν)
# Initial state vector
s0 = initState(oe_sat, oe_deb)

## Parameters
# Action space in m/s
A = [-10.0, 0.0, 10.0]
# MDP fields
proximityCost = -10
burnCost = -1
J₂ = 0
γ = 0.9
δt = 1
mdp = SATMDP(proximityCost, burnCost, Rₑ, μ, J₂, γ, A, δt)
# Forward Search parameters
depth = 5
steps = 5

## Run algorithim
sVec, rewards, actions, relDist = forwardSearch(mdp, depth, steps, s0)

## Display results
for i = 1 : steps
    println("Step ", i)
    println("Relative Distance: ", round(relDist[i]; digits=4), " km")
    println("Reward Received: ", round(rewards[i]; digits=4))
    println("Action Taken: ", actions[i])
    println(" ")
end
