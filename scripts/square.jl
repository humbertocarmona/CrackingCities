using OPC
using Geodesy
using Dates
using Statistics
using CSV, DataFrames
using Logging
using ArgParse

using Distributions
using Random

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--l1"
        arg_type = Float64
        default = 1000.0
        "--l2"
        arg_type = Float64
        default = 1000.0
        "--dl"
        arg_type = Float64
        default = 500.0
        "--nsamples", "-n"
        arg_type = Int
        default = 500
        "--nx"
        arg_type = Int
        default = 80
        "--runid"
        arg_type = String
        default = "1"
        "--prob"
        arg_type = Float64
        default = 0.6
        "--beta"
        arg_type = Float64
        default = -1.0
        "--logdir"
        arg_type = String
        default = "logs"
        "--resdir"
        arg_type = String
        default = "results"
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
# parsed_args = Dict{String,Any}(
#    "l1" => 1000.0,
#    "l2" => 1000.0,
#    "dl" => 1.0,
#    "nsamples" => 1,
#    "runid" => 1,
#    "prob" => 0.6,
#    "beta" => 0.06,
#    "nx" => 80,
#    "logdir" => "logs",
#    "resdir" => "results",
# )

l1 = parsed_args["l1"]
l2 = parsed_args["l2"]
dl = parsed_args["dl"]
ns = parsed_args["nsamples"]
runid = lpad(parsed_args["runid"],3,"0")
p = parsed_args["prob"]
β = parsed_args["beta"]
nx = parsed_args["nx"]
logdir = parsed_args["logdir"]
resdir = parsed_args["resdir"]

global dist = β
if β < 0
    # in case I want some other distribution I can read from a file
    df = CSV.read("../../data/boston-edges-4h.csv")|> DataFrame
    df[!,:tau] = df.tt ./ df.len
    n = size(df.tau,1)
    edges = collect(0:0.02:0.75)
    w = edges[2:end] - edges[1:end-1]
    dist = fit_mle(LogNormal, df.tau)
end

mkpath(logdir)
mkpath(resdir)

tinit = Dates.now()
tinits = Dates.format(tinit, "ddmmyy-HHhMM-SS")
io = open("$logdir/run$runid-$tinits.log", "w+")
logger = SimpleLogger(io, Logging.Info)
#global_logger(logger)
with_logger(logger) do
    @info("--------------------------
    $(Dates.format(tinit, "yy-mm-dd H:M"))
    runid = $runid, square lattice $nx x $nx, p = $p β = $β
    L = [$l1:$dl:$l2], nsamples = $ns
    --------------------------")
    flush(io)
    # if here, same network for each ℓ
    global weightmx = []
    for ℓ in collect(l1:dl:l2)
        # if here, one network for each ℓ
        nremoved = []
        measdist = []
        for sample = 1:ns
            # if here, one network for each od....
            seed = Dates.value(DateTime(Dates.now()))
            g, coords, weightmx, dd = OPC.buildSquareNetwork(
                nx,
                nx,
                p = p,
                dist = dist,
                seed = seed,
                od = false,
            )
            res_cellList = OPC.cellList(coords; cellWidth = 100.0)
            OD = OPC.odMatrix(
                ℓ,
                res_cellList;
                nDstOrg = 1,
                seed = seed,
                ns = 1,
                square = false,
                δ = 0.001,
            )

            (orig, dest) = OD[1]

            p1 = LLA(coords[orig][1], coords[orig][2], 0.0)
            p2 = LLA(coords[dest][1], coords[dest][2], 0.0)
            nrem, gr, rmmx = OPC.crackOptimalPaths(g, orig, dest, weightmx)
            # if sample == 1 && ℓ==1000.0
            #     fname = "results/l$(ℓ)-beta$(β)-p$(p).gpkg"
            #     OPC.writeShapeFile(g, coords, rmmx,fname)
            # end
            if mod(sample, 100) == 0
                @info("ℓ = $ℓ, sample = $sample, $(mean(nremoved))")
                flush(io)
            end
            push!(nremoved, nrem)
            push!(measdist, distance(p1, p2))
        end
        fn = "$resdir/$resdir-nr-$runid-l-$(Int(round(ℓ))).csv"
        df = DataFrame(ell = measdist, nr = nremoved)
        CSV.write(fn, df)
    end
    tend = Dates.now()
    dur = Dates.canonicalize(Dates.CompoundPeriod(tend - tinit))
    @info("----------- finished -----------------
    $(Dates.format(tend, "yy-mm-dd H:M"))
    took  $dur
    --------------------------------------------")
    flush(io)
end
close(io)

#vals = [i for i in findnz(weightmx)[3]]
#df = DataFrame(Dict("val"=>vals))
#CSV.write("vals.csv", df)
