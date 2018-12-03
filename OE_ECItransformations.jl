
"""
OE to ECI
"""
function OE2ECI(oe::OrbitalElementVec,μ=398600.4418)
    P = oe.a*(1-oe.e^2);                # Semi-Latus Rectum
    r_mag = P/(1+oe.e*cos(oe.ν));       # Distance from Earth to orbiting body

    n = sqrt(μ/oe.a^3)
    E = anom2E(oe.ν,oe.e)
    # R in perifocial coordinates [P Q W]'
    # r_peri = [r_mag*cos(oe.ν); r_mag*sin(oe.ν); 0];
    # v_peri = sqrt(μ/P)*[-sin(oe.ν); (oe.e+cos(oe.ν)); 0];
    r_peri = [oe.a*(cos(E) - oe.e); oe.a*sqrt(1 - oe.e^2)*sin(E);0];
    v_periComp = [-sin(E);sqrt(1 - oe.e^2)*cos(E);0];
    v_peri = (oe.a*n)/(1 - oe.e*cos(E))*v_periComp;
    if oe.i == 0 && oe.e != 0         # Equitorial and Elliptical
        R1 = 1;
        R2 = 1;
        R3 = rotz(oe.ω);
    elseif oe.e == 0 && oe.i != 0     # Circular and Inclined
        R1 = rotz(oe.Ω);
        R2 = rotx(oe.i);
        R3 = 1;
    elseif oe.i == 0 && oe.e == 0     # Equitorial and Circular
        R1 = 1;
        R2 = 1;
        R3 = 1;
    else                              # Not Circular or Inclined
        R1 = rotz(oe.Ω);
        R2 = rotx(oe.i);
        R3 = rotz(oe.ω);
    end
    R = R1*R2*R3;                     # Full rotation matrix
    r_eci = R*r_peri;
    v_eci = R*v_peri;
    return r_eci, v_eci
end



"""
ECI to OE
"""
function  ECI2OE(R,V,μ=398600.4418)
    # INPUTS
    #   R - I,J,K components of position
    #   V - I,J,K components of velocity
    # OUTPUTS
    #   oe - OrbitalElementVec
    # Function by
    #   Shawn Swist ~ 2018

    r = norm(R)
    v = norm(V)

    H = cross(R,V)
    h = norm(H)

    N = cross([0;0;1],H)
    n = norm(N)
    e_vec = 1/μ*((v^2-μ/r).*R-dot(R,V).*V)
    e = norm(e_vec)

    # Mechanical Energy to determine size
    ϵ = 0.5*v^2 - μ/r
    if e != 1
        a = -μ/(2*ϵ)
    else
        a = inf # Semi-major axis undefined for parabolas
    end

    # Orbital Inclination (always less than 180 deg)
    i = acos(H[3]/h)

    # Rignt Ascension of Ascending Node
    Ω = acos(N[1]/n)
    if N[2] < 0             # If Nⱼ is greater than 0 Om is less than 180
        Ω = 2π- Ω
    end

    # Argument of periapsis
    term = dot(N,e_vec)/(n*e)
    ϵ = 1e-10
    if abs(term) > 1 # checking precision of values
        if abs(term)-1 < ϵ
            if term < 0 term = -1 else term = 1 end
        end
    end
    ω = acos(term)
    if e_vec[3] < 0         # If e(k) is greater than 0 w is less than 180
        ω = 2π - ω;
    end

    # True anomaly
    term = dot(e_vec,R)/(e*r)
    ϵ = 1e-10
    if abs(term) > 1
        if abs(term)-1 < ϵ
            if term < 0 term = -1 else term = 1 end
        end
    end

    ν = acos(term);
    if dot(R,V) < 0         # If R dot V is greater than zero nu is less than 180
        ν = 2π - ν;
    end

    # Special Cases, use different angles
    if i == 0 && e != 0 # Elliptical equatorial
        # Provide the longitude of periapsis (PI = Om + w)
        ang = acos(e_vec[1]/e)
        if e_vec[2] < 0
            ang = 2π - ang;
        end
    elseif i != 0 && e == 0 # Circular inclined
        # Provide the argument of latitude (u = w + anom)
        ang = acos(dot(N,R)/(n*r))
        if r_eci[3] < 0
            ang = 2π - ang;
        end
    elseif i == 0 && e == 0 # Circular equatorial
        # Provide the true latitude (lambda = Om + w + anom)
        ang = acos(R[1]/r)
        if R[2] < 0
            ang = 2π - ang;
        end
    else
        # Default output for ang
        ang = NaN;
    end

    oe = OrbitalElementVec(a,e,i,Ω,ω,ν)
    return oe, ang
end
