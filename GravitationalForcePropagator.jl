"""
Gravitational Force propagation with J2
Run in ode solver
Inputs: time (for solver)
        X - [6x1] position and velocity state vector
        μ - Gravitational Parameter
        rₑ - Radius of center body
        J₂ - Zonal perturbation force
Output: Ẋ - Derivative of state vector
"""

function GravitationalForce(Ẋ,X,p,t)

    # Read μ, J₂, rₑ from parameters
    μ,rₑ,J₂ = p

    r = X[1:3]
    v = X[4:6]
    R = norm(r)
    r̂ = normalize(r)
    accel = -μ/R^2*r̂

    # Acceleration due to oblateness
    k̂ = [0 0 1]'
    z = r[3]
    accJ2 = - (μ*J₂*rₑ*rₑ/2)*((6*z/(R^5))*k̂ + ((3/(R^4)) - (15*z*z/(R^6)))*r̂)

    Ẋ[1:3] = v
    Ẋ[4:6] = accel+accJ2
end
