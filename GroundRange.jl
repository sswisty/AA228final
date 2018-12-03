"""
Find when a satellite is visible from a specific ground location
    do i want lat/long coords in degrees or radians???
"""


function GroundRange(Recef,ϕground,λground,rₑ=6378.137)

    # Convert latitude/longitude to radians
    ϕground = deg2rad.(ϕground)
    λground = deg2rad.(λground)
    Rgs = rₑ*MakeRgsUnitVec.(ϕground,λground)

    # Ground Station and Rotation Matrix
    Rxyzenu = ENUvec.(ϕground,λground)
    if length(Rgs) == length(Recef)
        Rsat = Rxyzenu.*(Recef.-Rgs)
    elseif length(Rgs) == 3
        Rsat = []
        for vect in Recef
            push!(Rsat,Rxyzenu*(vect-Rgs))
        end
    else
        error("Size of ground station does not match that of satellite. Must either be constant or have equivalent size.")
    end
    return Rsat
end

function SatAzAlt(Renu)
    re,rn,ru = [],[],[]
    for vect in Renu
        push!(re,vect[1])
        push!(rn,vect[2])
        push!(ru,vect[3])
    end
    Az = atan.(re,rn)
    Alt = atan.(ru,sqrt.(re.^2+rn.^2))
    # viz = findall(Alt.>0) # This gives locations where satellite is visible
    return Az, Alt
end

"""
Other functions that help vectorize code
"""

function MakeRgsUnitVec(ϕ,λ)
    r̂ = [cos(ϕ)*cos(λ);
         cos(ϕ)*sin(λ);
         sin(ϕ)]
    return r̂
end

function ENUvec(ϕ,λ)
    ê = [-sin(λ);cos(λ);0]
    n̂ = [-sin(ϕ)*cos(λ);-sin(ϕ)*sin(λ);cos(ϕ)]
    û = [cos(ϕ)*cos(λ);cos(ϕ)*sin(λ);sin(ϕ)]
    Rxyzenu = [ê n̂ û]'
    return Rxyzenu
end
