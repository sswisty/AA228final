"""
Basic time functions
    MJD2GMST
    UT12MJD

"""


function MJD2GMST(MJD)
    GMST = 280.4606.+360.9856473*(MJD.-51544.5);
    GMST = GMST*π/180;
    GMST = mod.(GMST,2π);
    return GMST
end



function UT12MJD( m, d, y, h, min, s )
#   Inputs
#       m/d/y - The date (month/day/year)
#       h:min:s - The time (24hour clock)
#
#   Returns
#       MJD - Modified Julian Date

    if m <= 2
        Y = y-1
        M = m+12
    else
        Y = y
        M = m
    end

    D = d + h/24 + min/1440 + s/86400
    B = Y/400 - Y/100 + Y/4
    MJD = 365*Y - 679004 + floor(B) + floor(30.6001*(M+1)) + D
    return MJD
end
