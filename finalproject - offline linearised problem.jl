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
using D3Trees

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
    px::Int64   #preferred position
    cr::Int64   #collision radius
end

function POMDPs.initialstate_distribution(p::SSAPOMDP)
    return Deterministic(p.start)
end

# Number of states, actions, observations
function POMDPs.n_states(p::SSAPOMDP)
    return p.nt*p.nt*p.nx*p.nx
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
        for ii = 1:p.nt
            for j = 1:p.nx
                for k = 1:p.nx
                    push!(x,[i, ii, j, k])
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
    return round(Int64,(s[1]-1)*p.nx*p.nx*p.nt + (s[2]-1)*p.nx*p.nx + (s[3]-1)*p.nx + s[4])
end
function POMDPs.actionindex(p::SSAPOMDP, a::Int64)
    return a
end
function POMDPs.obsindex(p::SSAPOMDP, o)
    return POMDPs.stateindex(p, o)
end

function POMDPs.discount(p::SSAPOMDP)
    return 0.999
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
    sp1 = [0;0;0;0]
    sp1[1] = s[1] + 1   # advance time
    sp1[2] = s[2]
    sp1[3] = s[3]       # debris position doesn't change
    sp1[4] = s[4] + p.actions[a]    # satellite position changes with manoeuvre
    sp2 = [0;0;0;0]
    sp2[1] = s[1] + 1
    sp2[2] = s[2]
    sp2[3] = s[3]
    sp2[4] = s[4]       # if manoeuvre fails, no change
    if sp1[4] > 10
        sp1[4] = 10     # edge of domain
    elseif sp1[4] < 1
        sp1[4] = 1      # edge of domain
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
    obs = [[0;0;0;0]]
    psat = [0.05;0.9]
    pdeb = [0.05;0.15;0.6]
    for i = -2:2
        if (sp[2] + i) < 1
            continue
        elseif (sp[2] + i) > p.nt
            continue
        else
            for j = -2:2
                if (sp[3] + j) < 1
                    continue
                elseif (sp[3] + j) > p.nx
                    continue
                else
                    for k = -1:1
                        if (sp[4] + k) < 1
                            continue
                        elseif (sp[4] + k) > p.nx
                            continue
                        else
                            push!(obs, sp + [0;i;j;k])
                            push!(probs, plist(p, i, j, k, sp, a))
                        end
                    end
                end
            end
        end
    end
    deleteat!(obs,1)
    deleteat!(probs,1)
    #return obs, probs
    return SparseCat(obs, probs)
end
function plist(p, t, d, s, sp, a)
    if a == 1
        p1 = phelp4(p, t, sp[2])
        p2 = phelp4(p, d, sp[3])
        p3 = phelp6(p, s, sp[4])
    else
        p1 = phelp1(p, t, sp[2])
        p2 = phelp1(p, d, sp[3])
        p3 = phelp3(p, s, sp[4])
    end
    return trunc(p1*p2*p3, digits = 10)
end
function phelp2(p, d, s)
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
function phelp1(p, t, s)
    if s >= 3 || s <= (p.nt-2)
        if abs(t) == 2
            pr = 0.05
        elseif abs(t) == 1
            pr = 0.15
        else
            pr = 0.6
        end
    elseif s == 2
        if t == -1
            pr = 0.2
        elseif t == 0
            pr = 0.6
        elseif t == 1
            pr = 0.15
        else
            pr = 0.05
        end
    elseif s == (p.nt-1)
        if t == 1
            pr = 0.2
        elseif t == 0
            pr = 0.6
        elseif t == -1
            pr = 0.15
        else
            pr = 0.05
        end
    else
        if t == 0
            pr = 0.75
        elseif abs(t) == 1
            pr = 0.2
        else
            pr = 0.05
        end
    end
    return pr
end
function phelp3(p, sat, s)
    if sat == 0
        p3 = 0.9
    else
        if s == 1 || s == p.nx
            p3 = 0.1
        else
            p3 = 0.05
        end
    end
    return p3
end
function phelp4(p, t, s)
    if s >= 3 || s <= (p.nt-2)
        if abs(t) == 2
            pr = 0.0
        elseif abs(t) == 1
            pr = 0.1
        else
            pr = 0.8
        end
    elseif s == 2
        if t == -1
            pr = 0.1
        elseif t == 0
            pr = 0.8
        elseif t == 1
            pr = 0.1
        else
            pr = 0.0
        end
    elseif s == (p.nt-1)
        if t == 1
            pr = 0.1
        elseif t == 0
            pr = 0.8
        elseif t == -1
            pr = 0.1
        else
            pr = 0.0
        end
    else
        if t == 0
            pr = 0.9
        elseif abs(t) == 1
            pr = 0.1
        else
            pr = 0.0
        end
    end
    return pr
end
function phelp6(p, t, s)
    if s >= 2 || s <= (p.nt-1)
        if abs(t) == 1
            pr = 0.05
        else
            pr = 0.9
        end
    elseif s == 1
        if t == -1
            pr = 0.0
        elseif t == 0
            pr = 0.95
        else
            pr = 0.05
        end
    else
        if t == -1
            pr = 0.05
        elseif t == 0
            pr = 0.95
        else
            pr = 0.0
        end
    end
    return pr
end

# Calculate reward
function POMDPs.reward(p::SSAPOMDP, s::Array{Int64,1}, a::Int64, sp::Array{Int64,1})
    r1 = 0.0
    distance = abs(sp[3] - sp[4])
    if sp[1] == sp[2]
        if distance <= p.cr
            r1 = r1 - 50        # cost of being within collision radius
        end
    end
    if a == 1
    elseif a <= 5
        r1 = r1 - 1*abs(p.actions[a]) # manoeuvre cost
    else
        r1 = r1 - 2*abs(p.actions[a])   # manoeuvre cost (more expensive to use N thruster)
    end
    if sp[4] == p.px
        r1 = r1 + 1.5         # reward for staying in nominal orbit
    end
    return r1
end

function init_pomdp()
    start = [1;20;5;3]
    nt = 50     # 15 timesteps
    nx = 10     # 10 position values
    actions = [0; -2; -1; 1; 2; -2; -1; 1; 2]   # thruster actions
    pt = 0.75    # probability of T thruster working
    pn = 0.95    # probability of N thruster working (more expensive but more reliable)
    px = 5      # nominal position
    cr = 1      # collision radius
    return SSAPOMDP(start, nt, nx, actions, pt, pn, px, cr)
end

function solve_pomdp()
    pomdp = init_pomdp()
    #solver = QMDPSolver(10000, 1e-8, false)
    #solver = FIBSolver()
    #solver = SARSOPSolver()
    #solver = DESPOTSolver(bounds=IndependentBounds(DefaultPolicyLB(RandomSolver()),3,check_terminal=true))
    #solver = POMCPSolver(max_depth = 40, c = 10.0, tree_queries = 5000)
    #solver = MCVISolver()
    solver = POMCPOWSolver(max_depth = 30, criterion = MaxUCB(10.0), tree_queries = 500)
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
    history = simulate(HistoryRecorder(max_steps = 50), pomdp, policy, belief_updater)
    for (s, b, a, o, r) in eachstep(history, "sbaor")
        # b = b.b
        println("State was $s,")
        #println("belief was $b,")
        println("action $a was taken,")
        println("the reward was $r,")
        println("and observation $o was received.\n")
    end
    #println("Discounted reward was $(discounted_reward(history)).")
    #a, info = action_info(policy, initialstate_distribution(pomdp), tree_in_info=true)
    #inchrome(D3Tree(info[:tree], init_expand=3))
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
