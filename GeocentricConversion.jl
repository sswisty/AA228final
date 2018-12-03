
function ECI2ECEF(R,V,θ)
    Recef = rotz.(-θ).*R
    Vecef = rotz.(-θ).*V

    return Recef,Vecef
end


function ECEF2GEO(R,rₑ=6378.137)
    #   Turn a position vector (in ECEF) into latitude and longitude
    #   Input
    #       R - 3x1 vector of position in ECEF frame
    #       rₑ - Radius of center body, defaults to earth
    #   OUTPUTS
    #       ϕ - latitude [degrees]
    #       λ - longitude [degrees]
    #       h - Altitude [km]
    r = norm.(R);
    rx,ry,rz,h = [],[],[],[]
    for vect in R
        push!(rx,vect[1])
        push!(ry,vect[2])
        push!(rz,vect[3])
        push!(h,norm(vect)-rₑ)
    end
    ϕ = asin.(rz./r)*180/π
    λ = atan.(ry,rx)*180/π
    return ϕ,λ,h
end

function LatLong2ECEF(ϕ,λ,rₑ=6378.137)
    # Note: Lat/Long are in degrees
    R = rₑ*[cosd(ϕ)*cosd(λ);cosd(ϕ)*sind(λ);sind(ϕ)]
    return R
end
