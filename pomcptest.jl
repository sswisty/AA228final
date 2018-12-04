using POMDPs
using POMCPOW
using POMDPModels
using POMDPSimulators
using POMDPPolicies

solver = POMCPOWSolver(criterion=MaxUCB(20.0))
pomdp = BabyPOMDP() # from POMDPModels
planner = POMDPs.solve(solver, pomdp)

hr = HistoryRecorder(max_steps=100)
hist = simulate(hr, pomdp, planner)
for (s, b, a, r, sp, o) in hist
    @show s, a, r, sp
end

rhist = simulate(hr, pomdp, RandomPolicy(pomdp))
println("""
    Cumulative Discounted Reward (for 1 simulation)
        Random: $(discounted_reward(rhist))
        POMCPOW: $(discounted_reward(hist))
    """)
