"""
Anomaly transformations: ν -> E -> M -> E -> ν
Convert true anomaly to eccentric anomaly to mean anomaly
Convert mean anomaly to eccentric anomaly to true anomaly
Mean to Eccentric anomaly has a tolerance, default set to 1e-8

"""

function anom2E(ν,e)
    E = acos((e + cos(ν))/(1 + e*cos(ν)));
    if ν > π
        E = 2π - E;
    end
    return E
end

function E2M(E,e)
    M = E - e*sin(E);
    return M
end

function M2E(M,e,tol=1e-8)
    if M == 0 || M == pi
        # Return the known solutions (trivial)
        E = M;
    else
        # Set up the problem based on an initial guess
        E0 = M;
        d = -(E0 - e*sin(E0) - M)/(1 - e*cos(E0));
        # Loop until the solution converges
        while abs(d) > tol
            E1 = E0 + d;
            d = -(E1 - e*sin(E1) - M)/(1 - e*cos(E1));
            E0 = E1;
        end
        E = E0;
    end
    return E
end

function E2anom(E,e)
    ν = acos((cos(E) - e)/(1 - e*cos(E)));
    if E > π
        ν = 2π - ν;
    end
    return ν
end
