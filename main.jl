include("OrbitMechFcns.jl")
include("MDPSAT_Helper.jl")

Rₑ = 6378.137
μ = 3.986e5
# Define the orbital elements of the satellite
a = Rₑ + 500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(0)
ω = 0
ν = deg2rad(0.001)
oe_sat = OrbitalElementVec(a,e,i,Ω,ω,ν)

# Define the orbital elements of the debris
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(0)
ω = 0
ν = deg2rad(0)
oe_deb = OrbitalElementVec(a,e,i,Ω,ω,ν)

A = [-10.0, 0.0, 10.0]
s0 = initState(oe_sat, oe_deb)
mdp = SATMDP(-10, -1, Rₑ, μ, 0, 0.9, A, 1)
depth = 10

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

# @time begin
# a, v = selectAction(s0, depth, mdp)
# end
# println(a, v)



action_taken = []
rsat = []
vsat = []
rdeb= []
si = []
rewards = []
push!(rsat,s0[1:3])
push!(vsat,s0[4:6])
push!(rdeb,s0[7:9])
push!(si,s0)

s = s0
stest = []
for i = 1:5
    push!(stest,s)
end

for i = 1:7

    a,v = selectAction(si[i],depth,mdp)
    sp = simulateDynamics(si[i], a, mdp)
    rew = rewardModel(si[i],a,mdp)
    push!(action_taken,a)
    push!(rsat,sp[1:3])
    push!(vsat,sp[4:6])
    push!(rdeb,sp[7:9])
    push!(rewards,rew)
    push!(si,sp)
end
