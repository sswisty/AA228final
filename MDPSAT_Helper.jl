struct SATMDP
    proximityCost::Float64
    burnCost::Float64
    Rₑ::Float64
    μ::Float64
    J₂::Float64
    discount::Float64
    A::Array{Float64,1}
    dt::Float64
end

function initState(oe_sat::OrbitalElementVec, oe_deb::OrbitalElementVec)
    Rs, Vs = OE2ECI(oe_sat)
    Rd, Vd = OE2ECI(oe_deb)
    s = [Rs; Vs; Rd; Vd]
    return s
end

function simulateDynamics(s, a, mdp::SATMDP)
    params = [mdp.μ, mdp.Rₑ, mdp.J₂]
    tspan = (0, mdp.dt)
    r = s[1:3]
    r̂ = normalize!(r)
    v = s[4:6]
    h = cross(r,v)
    n̂ = normalize!(h)
    t̂ = normalize!(cross(n̂,r̂))
    RotMat = [r̂ t̂ n̂]
    Δv = RotMat * [0; a / 1000; 0]
    IC1 = [s[1:3]; v + Δv]  # Position and Velocity of satellite
    IC2 = s[7:12] # Position and Velocity of debris
    prob1 = ODEProblem(GravitationalForce, IC1, tspan, params)
    prob2 = ODEProblem(GravitationalForce, IC2, tspan, params)
    sol1 = DifferentialEquations.solve(prob1, saveat = mdp.dt)
    sol2 = DifferentialEquations.solve(prob2, saveat = mdp.dt)
    # Build the next state vector
    sp = [sol1.u[end]; sol2.u[end]]
    return sp
end

function rewardModel(s, a, mdp::SATMDP)
    rBurn = mdp.burnCost * abs(a)
    relativeDistance = norm(s[1:3] - s[7:9])
    if relativeDistance > 30
        rProximity = 0.0
    elseif relativeDistance < 0.1
        rProximity = - 1e6
    else
        rProximity = mdp.proximityCost / relativeDistance
    end
    r = rProximity + rBurn
return r
end

function selectAction(s, d, mdp::SATMDP)
    if d == 0
        return 0.0, 0.0
    end
    astar = 0;
    vstar = -1e10
    for a in mdp.A
        v = rewardModel(s, a, mdp)
        sp = simulateDynamics(s, a, mdp)
        aprime, vprime = selectAction(sp, d - 1, mdp)
        v += mdp.discount * vprime
        if v > vstar
            vstar = v
            astar = a
        end
    end
return astar, vstar
end

function forwardSearch(mdp::SATMDP, depth, steps, s0)
    actions = []
    relDist = []
    stateVec = []
    rewards = []
    push!(relDist, norm(s0[1:3] - s0[7:9]))
    push!(stateVec, s0)
    for i = 1 : steps
        a, v = selectAction(stateVec[i], depth, mdp)
        sp = simulateDynamics(stateVec[i], a, mdp)
        r = rewardModel(stateVec[i], a, mdp)
        push!(actions, a)
        push!(relDist, norm(sp[1:3] - sp[7:9]))
        push!(rewards, r)
        push!(stateVec, sp)
    end
    return stateVec, rewards, actions, relDist
end
