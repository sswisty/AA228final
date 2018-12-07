using Printf
function write2matlab(file, data, m)
    datamat = zeros(length(data),m)
    for i = 1:length(data)
        for j = 1:m
            datamat[i,j] = data[i][j]
        end
    end
    open(file, "w") do io
        for i = 1:length(data)
            for j = 1:m
                @printf(io, "%s, ", datamat[i,j])
            end
            @printf(io, "\n")
        end
    end
    # return statemat
end
