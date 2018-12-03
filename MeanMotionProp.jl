"""
Mean motion satellite propagator

"""

function MeanMotionProp(oe::OrbitalElementVec,epoch,t_vec,flag="ECI",μ=398600.4418,rₑ=6378.137)

    n = sqrt(μ/oe.a^3)
    GMST = MJD2GMST(epoch.+t_vec)

    M0 = E2M(anom2E(oe.ν,oe.e),oe.e)
    M = M0.+n*t_vec*86400
    M = mod.(M,2π)

    Reci =  [zeros(3) for _ in 1:length(t_vec)]
    Veci =  [zeros(3) for _ in 1:length(t_vec)]
    Recef = [zeros(3) for _ in 1:length(t_vec)]
    Vecef = [zeros(3) for _ in 1:length(t_vec)]

    for i = 1:length(t_vec)
        E = M2E(M[i],oe.e,1e-10)    # Eccentric Anomaly [rad]
        ν = E2anom(E,oe.e)          # True Anomaly [rad]
        oeᵢ = OrbitalElementVec(oe.a,oe.e,oe.i,oe.Ω,oe.ω,ν)
        R,V = OE2ECI(oeᵢ)

        Reci[i] = R
        Veci[i] = V
        Recef[i] = rotz(-GMST[i])*R
        Vecef[i] = rotz(-GMST[i])*V
    end

    if flag == "ECEF"
        return Recef, Vecef
    elseif flag == "ECI"
        return Reci, Veci
    else
        error("Declare which frame! ECI or ECEF")
    end
end
