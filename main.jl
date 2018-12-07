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
ν = deg2rad(89.95)
oe_sat = OrbitalElementVec(a,e,i,Ω,ω,ν)

# Define the orbital elements of the debris
a = Rₑ+500          # 500 km altitude orbit
e = 0.0             #  circular orbit
i = deg2rad(90)     #  polar orbit
Ω = deg2rad(90)
ω = 0
ν = deg2rad(89.95)
oe_deb = OrbitalElementVec(a,e,i,Ω,ω,ν)

A = [-1.0, 0.0, 1.0]
s0 = initState(oe_sat, oe_deb)
mdp = SATMDP(-100, -2, Rₑ, μ, 0, 0.9, A, 1)
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

@time begin
a, v = selectAction(s0, depth, mdp)
end
println(a, v)
