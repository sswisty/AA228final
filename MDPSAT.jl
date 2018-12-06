using POMDPs
using MCTS
using Random
using POMDPModels
using POMDPSimulators
using POMDPPolicies

include("OrbitMechFcns.jl")

# Test to see if anything happens...
Rₑ = 6378.137
# Define the orbital elements of the satellite
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(90)
ω = 0
ν = 0
oe_sat = OrbitalElementVec(a,e,i,Ω,ω,ν)
Rs,Vs = OE2ECI(oe_sat)

# Define the orbital elements of the debris
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(0)
ω = 0
ν = 0
oe_deb = OrbitalElementVec(a,e,i,Ω,ω,ν)
Rd,Vd = OE2ECI(oe_deb)

initialState = [Rs; Vs; Rd; Vd]

struct SatMDP <: MDP{Array{Float64,1}, Float64}
    r_close::Float64 # = -100      # 1/r^2 reward coefficient
    r_burn::Float64 # = -10        # burn reward, multiple for multiple directions?
    γ::Float64  # .9              # Discount factor, necessary?
    μ::Float64  # 398600.441      #
    Rₑ::Float64  # 6378.137       # [km] Earth radius
    J₂::Float64  # 1.0826267e-3   # Oblatness
    dt::Float64
end

POMDPs.actions(::SatMDP) = [-2.0 -1.0 0.0 1.0 2.0]

POMDPs.discount(::SatMDP) = 1

function POMDPs.generate_sr(P::SatMDP, s::Array{Float64,1}, a::Float64, rng::AbstractRNG) # s,a are arrays...
    # Define the constants for the propagator
    tspan = (0, P.dt)
    params = [P.μ, P.Rₑ, P.J₂]

    # Set up as 2 solutions - sat and debris
    IC1 = s[1:6]  # Position and Velocity of satellite
    IC2 = s[7:12] # Position and Velocity of debris
    prob1 = ODEProblem(GravitationalForce, IC1, tspan, params)
    prob2 = ODEProblem(GravitationalForce, IC2, tspan, params)
    sol1 = DifferentialEquations.solve(prob1, saveat = P.dt) # hope this is right ...
    sol2 = DifferentialEquations.solve(prob2, saveat = P.dt)

    # Build the next state vector
    sp = [sol1.u[end]; sol2.u[end]] # transpose keeps sp a row vector...

    # Define the reward
    r = norm(sp[1:3] - sp[7:9])
    Rclose = P.r_close*1/r^2 # maybe make this zero if below a threshold ...
    Rburn = P.r_burn*abs(a)   # if just tangential burns are possible, may redefine later

    R = Rclose + Rburn

    return sp,R
end

@requirements_info MCTSSolver() SatMDP(1,2,3,4,7,5,6) 1

# POMDPs.initialstate_distribution(p::SatMDP) = s
# satPOMDP = SatPOMDP(-100,-10,.9,398600.441,6378.137,1.0826267e-3,1)

mdp = SatMDP(-100, -10, 0.9, 398600.441, 6378.137, 1.0826267e-3, 1) # initializes the MDP
solver = MCTS.DPWSolver(n_iterations=50, depth=20, exploration_constant=5.0) # initializes the Solver type
policy = POMDPs.solve(solver, mdp) # initializes the planner

# for (s,a,r) in stepthrough(mdp, planner, "s,a,r", max_steps = 10)
#     @show s, a, r
# end

# state = initialstate(mdp, Random.MersenneTwister(4))
# a = action(policy, state)

# s = initialstate(mdp, Random.seed!(1))
# simulate(simulator::Simulator, problem::MDP{S,A}, policy::Policy, initial_state::S)
rsat = []
rdeb =[]
for (s, a, r) in stepthrough(mdp, policy, initialState, "s,a,r", max_steps=5)
    @show a, r, s
    push!(rsat,s[1:3])
    push!(rdeb,s[7:9])
end