using POMDPs
using Random # for AbstractRNG
using POMDPModelTools # for Deterministic

struct BabyPOMDP <: POMDP{Bool, Bool, Bool}
    r_feed::Float64
    r_hungry::Float64
    p_become_hungry::Float64
    p_cry_when_hungry::Float64
    p_cry_when_not_hungry::Float64
    discount::Float64
end

BabyPOMDP() = BabyPOMDP(-5., -10., 0.1, 0.8, 0.1, 0.9);

function POMDPs.generate_s(p::BabyPOMDP, s::Bool, a::Bool, rng::AbstractRNG)
    if a # feed
        return false
    elseif s # hungry
        return true
    else # not hungry
        return rand(rng) < p.p_become_hungry
    end
end
function POMDPs.generate_o(p::BabyPOMDP, s::Bool, a::Bool, sp::Bool, rng::AbstractRNG)
    if sp # hungry
        return rand(rng) < p.p_cry_when_hungry
    else # not hungry
        return rand(rng) < p.p_cry_when_not_hungry
    end
end
POMDPs.reward(p::BabyPOMDP, s::Bool, a::Bool) = s*p.r_hungry + a*p.r_feed
POMDPs.initialstate_distribution(m::BabyPOMDP) = Deterministic(false)

using POMDPSimulators
using POMDPPolicies

m = BabyPOMDP()

# policy that maps every input to a feed (true) action
policy = FunctionPolicy(o->true)

for (s, a, r) in stepthrough(m, policy, "s,a,r", max_steps=10)
    @show s
    @show a
    @show r
    println()
end
