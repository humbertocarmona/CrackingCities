function findFilesMaching(pathern::Regex, path::String)
    f = readdir(path)

    idx = findall(x->occursin(pathern,x), f)
    fs = f[idx]
    files = []
    for x in fs
        push!(files,"$path/$x")
    end
    files
end

function hist(data::Array{Int64,1}, edges::Array{Float64,1}; norm::Bool=false)
    d = convert.(Float64, data)
    h, c = hist(d, edges, norm=norm)
    return h,c
end

function hist(data::Array{Float64,1}, edges::Array{Float64,1}; norm::Bool=false)
    n = size(data,1)
    w = edges[2:end] - edges[1:end-1]
    c = 0.5*(edges[1:end-1]+edges[2:end])
    h = zeros(size(edges,1)-1)
    for val in data
        idx = findfirst(x->x>val, edges)
        if idx != nothing && idx > 1
            h[idx-1] += 1.0
        end
    end
    if norm
        h = h./w/n
    end
    return h, c
end

# using CSV, DataFrames
#
# files = findFilesMaching(r"square-nr-.+-l-.+.csv", "results")
#
#
# for f in files
#     df = CSV.read(f)|> DataFrame
#     unique!(df)
#     fl = replace(f, r"square" => s"c-square")
#     println("$fl")
#     println(first(df,5))
#     CSV.write(fl, df)
# end
