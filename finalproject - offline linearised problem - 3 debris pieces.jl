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
using BasicPOMCP# works, but worse result
using POMDPSolve# can't simulate?
using POMCPOW   # works, but worse result
using AEMS      # too slow (uses FIB Solve)

# state: [time, satellite position, debris 1 pos, debris 2 pos, debris 3 pos]
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
    ct::Array{Int64,1}   #collision times
    px::Int64   #preferred position
    cr::Array{Int64,1}   #collision radii
    nd::Int64   #number of debris pieces
end

function POMDPs.initialstate_distribution(p::SSAPOMDP)
    return Deterministic(p.start)
end

# Number of states, actions, observations
function POMDPs.n_states(p::SSAPOMDP)
    return p.nt*p.nx^(p.nd+1)
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
                for l = 1:p.nx
                    for m = 1:p.nx
                        push!(x,[i, j, k, l, m])
                    end
                end
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
    index = 0
    for i = 1:(p.nd+1)
        index = index + (s[i]-1)*p.nx^(5-i)
    end
    index = index + s[5]
    return round(Int64,index)
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
    sp1 = [0;0;0;0;0]
    sp2 = [0;0;0;0;0]
    sp1[1] = s[1] + 1   # advance time
    sp2[1] = s[1] + 1
    sp1[2] = s[2] + p.actions[a] # satellite position changes with manoeuvre
    sp2[2] = s[2]                # if manoeuvre fails, no change
    for i = 3:(p.nd+2)
        sp1[i] = s[i]   # debris position doesn't change
        sp2[i] = s[i]
    end
    if sp1[2] > 10
        sp1[2] = 10     # edge of domain
    elseif sp1[2] < 1
        sp1[2] = 1      # edge of domain
    end
    if a <= 5
        p = [p.pt; 1-p.pt]  # probabilities for first thruster
    else
        p = [p.pn; 1-p.pn]  # probabilities for second thruster
    end
    return SparseCat([sp1, sp2], p)
end

# Generate list of possible observations and associated probabilities
# Debris position has a +/-2 uncertainty range, p_correct = 0.6
# Satellite position has a +/-1 uncertainty range, p_correct = 0.9
function POMDPs.observation(p::SSAPOMDP, a::Int64, sp::Array{Int64,1})
    probs = [0.0]
    obs = [[0;0;0;0;0]]
    psat = [0.05;0.9]
    pdeb = [0.05;0.15;0.6]
    for i = -1:1
        if (sp[2] + i) < 1
            continue
        elseif (sp[2] + i) > p.nx
            continue
        else
            for j = -2:2
                if (sp[3] + j) < 1
                    continue
                elseif (sp[3] + j) > p.nx
                    continue
                else
                    for k = -2:2
                        if (sp[4] + k) < 1
                            continue
                        elseif (sp[4] + k) > p.nx
                            continue
                        else
                            for l = -2:2
                                if (sp[5] + l) < 1
                                    continue
                                elseif (sp[5] + l) > p.nx
                                    continue
                                else
                                    push!(obs, sp + [0;i;j;k;l])
                                    push!(probs, plist(p, i, j, k, l, sp))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    pop!(obs)
    pop!(probs)
    #return obs, probs
    return SparseCat(obs, probs)
end
function plist(p, s, d1, d2, d3, sp)
    if s == 0
        p1 = 0.9
    else
        if sp[2] == 1 || sp[2] == p.nx
            p1 = 0.1
        else
            p1 = 0.05
        end
    end
    p2 = phelp(p, d1, sp)
    p3 = phelp(p, d2, sp)
    p4 = phelp(p, d3, sp)
    return p1*p2*p3*p4
end
function phelp(p, d, s)
    if d == 0
        pr = 0.6
    elseif d == -1 && s == 2
        pr = 0.2
    elseif d == 1 && s == (p.nx-1)
        pr = 0.2
    elseif abs(d) == 1
        pr = 0.15
    elseif d == 0 && s == 1
        pr = 0.75
    elseif d == 1 && s == 1
        pr = 0.2
    elseif d == 0 && s == p.nx
        pr = 0.75
    elseif d == -1 && s == p.nx
        pr = 0.2
    else
        pr = 0.05
    end
    return pr
end

# Calculate reward
function POMDPs.reward(p::SSAPOMDP, s::Array{Int64,1}, a::Int64, sp::Array{Int64,1})
    r = 0.0
    time = sp[1]
    for i = 3:(p.nd+2)
        distance = abs(sp[2] - sp[i])
        if time == p.ct[i-2]
            if distance <= p.cr[i-2]
                r = r - 20  # cost of being within collision radius
            end
        end
    end
    if a == 1
    elseif a <= 5
        r = r - abs(p.actions[a]) # manoeuvre cost
    else
        r = r - 2*abs(p.actions[a])   # manoeuvre cost (more expensive to use N thruster)
    end
    if sp[2] == p.px
        r = r + 1         # reward for staying in nominal orbit
    end
    return r
end

function init_pomdp()
    start = [1;3;5;5;3]
    nt = 20     # 15 timesteps
    nx = 10     # 10 position values
    actions = [0; -2; -1; 1; 2; -2; -1; 1; 2]   # thruster actions
    pt = 0.7    # probability of T thruster working
    pn = 0.9    # probability of N thruster working (more expensive but more reliable)
    ct = [10; 12; 15]     # collision timestep
    px = 5      # nominal position
    cr = [1; 1; 1]      # collision radius
    nd = 3
    return SSAPOMDP(start, nt, nx, actions, pt, pn, ct, px, cr,nd)
end

function solve_pomdp()
    pomdp = init_pomdp()
    solver = QMDPSolver(10000, 1e-8, false)
    #solver = FIBSolver()
    #solver = SARSOPSolver()
    #solver = DESPOTSolver(bounds=IndependentBounds(DefaultPolicyLB(RandomSolver()),3,check_terminal=true))
    #solver = POMCPSolver()
    #solver = MCVISolver()
    #solver = POMCPOWSolver()
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
    history = simulate(HistoryRecorder(max_steps = 21), pomdp, policy, belief_updater)
    for (s, b, a, o, r) in eachstep(history, "sbaor")
        # b = b.b
        println("State was $s,")
        #println("belief was $b,")
        println("action $a was taken,")
        println("the reward was $r,")
        println("and observation $o was received.\n")
    end
    #println("Discounted reward was $(discounted_reward(history)).")
    return discounted_reward(history)
end

function sim_lots(runs)
    r = 0
    for i = 1:runs
        r = r + sim_pomdp()
    end
    r = r/runs
    return r
end

#rand_policy = RandomPolicy(pomdp)
