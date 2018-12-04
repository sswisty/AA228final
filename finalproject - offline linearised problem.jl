using POMDPs
using POMDPModels
using POMDPSimulators
using POMDPPolicies
using POMDPModelTools
using POMDPFiles
using Random
using LinearAlgebra
using BeliefUpdaters
using FIB       # too slow
using SARSOP    # crashes
using MCVI      # bad documentation
using QMDP      # works, but variable result
using ARDESPOT  # works, but weird policy sometimes
using BasicPOMCP# works, but not as good?
using POMDPSolve# can't simulate?
using POMCPOW   # can't simulate?
using AEMS      # too slow (uses FIB Solve)

# state: [time, debris position, satellite position]
# time is in discrete steps and advances at every state transition (e.g t = 1->2->3->4...)
# positions are in discrete steps (e.g. [1,2,3,4,5,6,7,8,9,10])
# I guess you could say that if the satellite/debris are close enough, it's basically like linearizing the dynamics?

# define structure to contain necessary info
struct SSAPOMDP <: POMDP{Array{Int64,1}, Int64, Array{Int64,1}}
    start::Array{Int64,1}   # initial state
    nt::Int64       # number of times
    nx::Int64       # number of distances
    actions::Array{Int64,1} # list of actions
    pt::Float64   #probability of T thruster working
    pn::Float64   #probability of N thruster working
    ct::Int64   #collision time
    px::Int64   #preferred position
    cr::Int64   #collision radius
end

function POMDPs.initialstate_distribution(p::SSAPOMDP)
    return Deterministic(p.start)
end

# Number of states, actions, observations
function POMDPs.n_states(p::SSAPOMDP)
    return p.nt*p.nx*p.nx
end
function POMDPs.n_actions(p::SSAPOMDP)
    return length(p.actions)
end
function POMDPs.n_observations(p::SSAPOMDP)
    return POMDPs.n_states(p)
end

# List of states, actions, observations
function POMDPs.states(p::SSAPOMDP)
    x = []
    for i = 1:p.nt
        for j = 1:p.nx
            for k = 1:p.nx
                push!(x,[i, j, k])
            end
        end
    end
    return x
end
function POMDPs.actions(p::SSAPOMDP)
    return [i for i in 1:length(p.actions)]
end
function POMDPs.observations(p::SSAPOMDP)
    return POMDPs.states(p)
end

# Convert state, action, observation to integer index
function POMDPs.stateindex(p::SSAPOMDP, s)
    return round(Int64,(s[1]-1)*p.nx*p.nx + (s[2]-1)*p.nx + s[3])
end
function POMDPs.actionindex(p::SSAPOMDP, a::Int64)
    return a
end
function POMDPs.obsindex(p::SSAPOMDP, o)
    return POMDPs.stateindex(p, o)
end

function POMDPs.discount(p::SSAPOMDP)
    return 0.95
end

# Is the state terminal (i.e. is t = t_end)
function POMDPs.isterminal(p::SSAPOMDP, s::Array{Int64,1})
    if s[1] == p.nt
        return true
    else
        return false
    end
end

function POMDPs.transition(p::SSAPOMDP, s::Array{Int64,1}, a::Int64)
    sp1 = [0;0;0]
    sp1[1] = s[1] + 1   # advance time
    sp1[2] = s[2]       # debris position doesn't change
    sp1[3] = s[3] + p.actions[a]    # satellite position changes with manoeuvre
    sp2 = [0;0;0]
    sp2[1] = s[1] + 1
    sp2[2] = s[2]
    sp2[3] = s[3]       # if manoeuvre fails, no change
    if sp1[3] > 10
        sp1[3] = 10     # edge of domain
    elseif sp1[3] < 1
        sp1[3] = 1      # edge of domain
    end
    if a <= 5
        p = [p.pt; 1-p.pt]  # probability for first thruster
    else
        p = [p.pn; 1-p.pn]  # probability for second thruster
    end
    return SparseCat([sp1, sp2], p)
end

# Generate list of possible observations and associated probabilities
# Debris position has a +/-2 uncertainty range, p_correct = 0.6
# Satellite position has a +/-1 uncertainty range, p_correct = 0.9
function POMDPs.observation(p::SSAPOMDP, a::Int64, sp::Array{Int64,1})
    probs = []
    states = []
    if sp[2] >= 3 || sp[2] <= 8
        probs, states = probs3(0.05, [sp[1]; sp[2]-2; sp[3]], probs, states)
        probs, states = probs3(0.15, [sp[1]; sp[2]-1; sp[3]], probs, states)
        probs, states = probs3(0.6, [sp[1]; sp[2]; sp[3]], probs, states)
        probs, states = probs3(0.15, [sp[1]; sp[2]+1; sp[3]], probs, states)
        probs, states = probs3(0.05, [sp[1]; sp[2]+2; sp[3]], probs, states)
    elseif sp[2] == 2
        probs, states = probs3(0.2, [sp[1]; sp[2]-1; sp[3]], probs, states)
        probs, states = probs3(0.6, [sp[1]; sp[2]; sp[3]], probs, states)
        probs, states = probs3(0.15, [sp[1]; sp[2]+1; sp[3]], probs, states)
        probs, states = probs3(0.05, [sp[1]; sp[2]+2; sp[3]], probs, states)
    elseif sp[2] == 9
        probs, states = probs3(0.05, [sp[1]; sp[2]-2; sp[3]], probs, states)
        probs, states = probs3(0.15, [sp[1]; sp[2]-1; sp[3]], probs, states)
        probs, states = probs3(0.6, [sp[1]; sp[2]; sp[3]], probs, states)
        probs, states = probs3(0.2, [sp[1]; sp[2]+1; sp[3]], probs, states)
    elseif sp[2] == 1
        probs, states = probs3(0.75, [sp[1]; sp[2]; sp[3]], probs, states)
        probs, states = probs3(0.2, [sp[1]; sp[2]+1; sp[3]], probs, states)
        probs, states = probs3(0.05, [sp[1]; sp[2]+2; sp[3]], probs, states)
    else
        probs, states = probs3(0.05, [sp[1]; sp[2]-2; sp[3]], probs, states)
        probs, states = probs3(0.2, [sp[1]; sp[2]-1; sp[3]], probs, states)
        probs, states = probs3(0.75, [sp[1]; sp[2]; sp[3]], probs, states)
    end
    probs = trunc.(probs, digits = 5)
    return SparseCat(states, probs)
end
function probs3(p::Float64, sp::Array{Int64,1}, probs, states)
    if sp[3] >= 2 || sp[3] <= 9
        push!(probs, 0.05*p)
        push!(probs, 0.9*p)
        push!(probs, 0.05*p)
        push!(states, [sp[1];sp[2];sp[3]-1])
        push!(states, [sp[1];sp[2];sp[3]])
        push!(states, [sp[1];sp[2];sp[3]+1])
    elseif sp[3] == 1
        push!(probs, 0.9*p)
        push!(probs, 0.1*p)
        push!(states, [sp[1];sp[2];sp[3]])
        push!(states, [sp[1];sp[2];sp[3]+1])
    else
        push!(probs, 0.1*p)
        push!(probs, 0.9*p)
        push!(states, [sp[1];sp[2];sp[3]-1])
        push!(states, [sp[1];sp[2];sp[3]])
    end
    return probs, states
end

# Calculate reward
function POMDPs.reward(p::SSAPOMDP, s::Array{Int64,1}, a::Int64, sp::Array{Int64,1})
    r1 = 0
    distance = abs(sp[2] - sp[3])
    time = sp[1]
    if time == p.ct
        if distance <= p.cr
            r1 = r1 - 20        # cost of being within collision radius
        end
    end
    if a == 1
    elseif a <= 5
        r1 = r1 - abs(p.actions[a]) # manoeuvre cost
    else
        r1 = r1 - 2*abs(p.actions[a])   # manoeuvre cost (more expensive to use N thruster)
    end
    if sp[3] == p.px
        r1 = r1 + 1         # reward for staying in nominal orbit
    end
    return r1
end

function init_pomdp()
    start = [1;5;3]
    nt = 15     # 15 timesteps
    nx = 10     # 10 position values
    actions = [0; -2; -1; 1; 2; -2; -1; 1; 2]   # thruster actions
    pt = 0.7    # probability of T thruster working
    pn = 0.9    # probability of N thruster working (more expensive but more reliable)
    ct = 10     # collision timestep
    px = 5      # nominal position
    cr = 1      # collision radius
    return SSAPOMDP(start, nt, nx, actions, pt, pn, ct, px, cr)
end

function solve_pomdp()
    pomdp = init_pomdp()
    solver = QMDPSolver(10000, 1e-8, false)
    #solver = FIBSolver()
    #solver = SARSOPSolver()
    #solver = DESPOTSolver(bounds=IndependentBounds(DefaultPolicyLB(RandomSolver()),3,check_terminal=true))
    #solver = POMCPSolver()
    #solver = MCVISolver()
    #solver = POMCPOWSolver(criterion=MaxUCB(20.0))
    #solver = AEMSSolver()
    policy = solve(solver, pomdp)
    #fout = open("SSA.pomdp", "w")
    #write(fout, pomdp)
    #close(fout)
    #pomdpf = POMDPSolveFile("SSA.pomdp")
    #solver = POMDPSolveSolver()
    #policy = POMDPSolvePolicy("SSA.policy")
    #policy = POMDPSolve.solve(solver, pomdpf, policy)

    return pomdp, policy
end

function sim_pomdp()
    pomdp, policy = solve_pomdp()
    belief_updater = updater(policy)
    history = simulate(HistoryRecorder(max_steps = 16), pomdp, policy, belief_updater)
    for (s, b, a, o, r) in eachstep(history, "sbaor")
        # b = b.b
        println("State was $s,")
        #println("belief was $b,")
        println("action $a was taken,")
        println("the reward was $r,")
        println("and observation $o was received.\n")
    end
    println("Discounted reward was $(discounted_reward(history)).")
    #return history
end

#rand_policy = RandomPolicy(pomdp)
