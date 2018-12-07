include("OrbitMechFcns.jl")
include("MDPSAT_Helper.jl")
include("export2matlab.jl")

## Constants
Rₑ = 6378.137
μ = 3.986e5

## Orbit Initialization
h = 500 # altitude
# Orbital elements of the satellite
a = Rₑ + h; e = 0.0; i = deg2rad(22); Ω = deg2rad(0); ω = 0; ν = deg2rad(0)
oe_sat = OrbitalElementVec(a, e, i, Ω, ω, ν)
# Orbital elements of the debris
a = Rₑ + h; e = 0.0; i = deg2rad(22); Ω = deg2rad(0); ω = 0; ν = deg2rad(0.001)
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
depth = 10
steps = 20

## Run algorithim
@time begin
sVec, rewards, actions, relDist = forwardSearch(mdp, depth, steps, s0)
end

## Display results
for i = 1 : steps
    println("Step ", i)
    println("Relative Distance: ", round(relDist[i]; digits=4), " km")
    println("Reward Received: ", round(rewards[i]; digits=4))
    println("Action Taken: ", actions[i])
    println(" ")
end

# Define the nominal orbits for plotting
ICs = s0[1:6]
ICd = s0[7:12]
tspan = (0.0,95*60)
params = [μ,Rₑ,J₂]
sat_ode = ODEProblem(GravitationalForce, ICs, tspan, params)
deb_ode = ODEProblem(GravitationalForce, ICd, tspan, params)
sol1 = DifferentialEquations.solve(sat_ode, saveat = 2)
sol2 = DifferentialEquations.solve(deb_ode, saveat = 2)

# Take the final state and make the new orbit
ICf = sVec[end][1:6]
sat_ode_f = ODEProblem(GravitationalForce, ICf, tspan, params)
sol3 = DifferentialEquations.solve(sat_ode_f, saveat = 2)

write2matlab("sat_nominal.csv",sol1.u,6)
write2matlab("deb_nominal.csv",sol2.u,6)
write2matlab("sat_end.csv",sol3.u,6)
write2matlab("states.csv",sVec,12)
