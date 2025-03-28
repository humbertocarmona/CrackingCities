using CSV
using DataFrames

edf = CSV.read("data/boston-edges-4h.csv")
ndf = CSV.read("data/boston-nodes.csv")

nn = size(ndf,1)
nid = Dict()
for i in 1:nn
    id = ndf.id[i]
    nid[id] = i
end

newids = [nid[ndf.id[i]] for i in 1:nn]
oldids = [ndf.id[i] for i in 1:nn]
lat = [ndf.lat[i] for i in 1:nn]
lon = [ndf.lon[i] for i in 1:nn]
df = DataFrame(i=newids, n=oldids, lat=lat, lon=lon)

# CSV.write("data/boston-nodes.csv", df)


ne = size(edf,1)
src = [nid[edf.source[i]] for i in 1:ne]
dst = [nid[edf.target[i]] for i in 1:ne]
dist = [edf.dis[i] for i in 1:ne]
wgt = [edf.tt[i] for i in 1:ne]

df = DataFrame(src=src, dst=dst, len=dist, tt=wgt)
sort!(df, (:src, :dst), rev=(false, false))
# CSV.write("data/boston-edges-4h.csv", df)



#  using JLD, HDF5
# function getWeights(fname = "data/fortaleza_traveltime.jld")
#     traveltime = load(fname, "traveltime")
#     N = nv(g)
#     weightmx = spzeros(N,N)
#     for e in edges(g)
#         i, j = e.src, e.dst
#         weightmx[i,j] = 0.06
#         τ = traveltime[(i,j)]/distmx[i,j]
#         if τ > 0.06
#             weightmx[i,j] = τ
#         end
#     end
#     return weightmx
# end
#
#
# close(io)
