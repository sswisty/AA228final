

function EarthGroundPlot()
    topo = CSV.read("/Users/SSwisty/Documents/Julia/OrbitalMechanics/topo.csv",header=false)
    # topo = CSV.read("D:/Documents/Julia/OrbitalMechanicsFcns/topo.csv",header=false)
    n,m = size(topo)
    data =zeros(n,m)
    for i = 1:n
        for j = 1:m
            data[i,j] = topo[i,j]
        end
    end
    mapvals = [data[:, 181:360] data[:, 1:180]]
    plt = contour(-180:179, -90:89, mapvals; levels=0:1,legend=false)
    return plt
end
