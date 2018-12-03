"""
AA 228 Final Project
Justin, Andrew, Shawn
Stanford University, Fall 2018

NOTES
I have been setting this up like we talked about and using
https://github.com/JuliaPOMDP/POMDPExamples.jl/blob/master/notebooks/Defining-a-POMDP-with-the-Generative-Interface.ipynb
as a reference.
"""

# Packages used
using POMDPs
using POMDPModels               # necessary?
using POMDPSimulators           # necessary?
using DifferentialEquations
include("OrbitMechFcns.jl")

# Defining the POMDP
struct SatPOMDP <: POMDP{float64,float64,float64} # The three inputs are for state, action, observation, all floats
    r_close::Float64 # = -100      # 1/r^2 reward coefficient
    r_burn::Float64 # = -10        # burn reward, multiple for multiple directions?
    γ::Float64 = # .9              # Discount factor, necessary?
    μ::Float64 = # 398600.441      #
    Rₑ::Float64 = # 6378.137       # [km] Earth radius
    J₂::Float64 = # 1.0826267e-3   # Oblatness
    dt::Float64
end

satPOMDP = SatPOMDP(-100,-10,.9,398600.441,6378.137,1.0826267e-3,1)

function POMDPs.generate_sor(P::SatPOMDP, s, a) # s,a are arrays...
    # Define the constants for the propagator
    tspan = (0,P.dt)
    params = [P.μ,P.Rₑ,P.J₂]

    # Set up as 2 solutions - sat and debris
    IC1 = s[1:6]  # Position and Velocity of satellite
    IC2 = s[7:12] # Position and Velocity of debris
    prob1 = ODEProblem(GravitationalForce,IC1,tspan,params)
    prob2 = ODEProblem(GravitationalForce,IC1,tspan,params)
    sol1 = DifferentialEquations.solve(prob1,saveat = dt) # hope this is right ...
    sol2 = DifferentialEquations.solve(prob2,saveat = dt)

    # Build the next state vector
    sp = [sol1.u[end] sol2.u[end]] # make sure this pulls out all the right values

    # Define the reward
    r = abs(sp[1:3]-sp[7:9]) # probably an array size issue here
    Rclose = P.r_close*1/r^2
    Rburn = P.r_burn*abs(a)   # if just tangential burns are possible, may redefine later

    R = Rclose + Rburn

    # OBSERVATION ...
    o = sp      # for testing...


    return sp,R,o
end
