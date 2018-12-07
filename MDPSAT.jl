using POMDPs
using MCTS
using Random
using POMDPModels
using POMDPSimulators
using POMDPPolicies

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

    # Rotate the action (Δv) from RTN frame into ECI frame for propagator




    return sp,R
end

# @requirements_info MCTSSolver() SatMDP(1,2,3,4,7,5,6) 1


mdp = SatMDP(-1000, -10, 0.9, 398600.441, 6378.137, 1.0826267e-3, 1)
solver = MCTS.DPWSolver(n_iterations=50, depth=20, exploration_constant=0.0)
policy = POMDPs.solve(solver, mdp) # initializes the planner

rsat = []
rdeb =[]
approach = []
action_taken = []
reward_recieved = []
for (s, a, r) in stepthrough(mdp, policy, initialState, "s,a,r", max_steps=100)
    # @show a, r, s
    push!(rsat,s[1:3])
    push!(rdeb,s[7:9])
    push!(approach,norm(s[1:3]-s[7:9]))
    push!(action_taken,a)
    push!(reward_recieved,r)
end

plot(approach)
xsat = [];ysat = [];zsat = []
for vect in rsat
    push!(xsat,vect[1])
    push!(ysat,vect[2])
    push!(zsat,vect[3])
end
xdeb = [];ydeb = [];zdeb = []
for vect in rdeb
    push!(xdeb,vect[1])
    push!(ydeb,vect[2])
    push!(zdeb,vect[3])
end

plot(xsat,ysat,zsat)
plot!(xdeb,ydeb,zdeb)
