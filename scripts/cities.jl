using OPC
using Geodesy
# using LightGraphs
# using SparseArrays
using Dates
using Statistics
using CSV, DataFrames
using Logging
using ArgParse

# googlekey = "AIzaSyApQzC_OLdxiITS7ynh_XsWZZOU8XOKQHs"

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
            default = 1000.0
        "--nsamples"
            arg_type = Int
            default = 500
        "--runid"
            arg_type = String
            default = "1"
        "--efile"
            arg_type = String
            default = ""
        "--nfile"
            arg_type = String
            default = ""
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
#    "nsamples" => 1000,
#    "runid" => "1",
#    "efile"=>"../data/boston-edges-4h.csv",
#    "nfile"=>"../data/boston-nodes.csv",
#    "logdir" => "logs",
#    "resdir" => "results",
# )

l1 = parsed_args["l1"]
l2 = parsed_args["l2"]
dl = parsed_args["dl"]
ns = parsed_args["nsamples"]
runid = lpad(parsed_args["runid"],3,"0")
efile = parsed_args["efile"]
nfile = parsed_args["nfile"]
logdir = parsed_args["logdir"]
resdir = parsed_args["resdir"]
city = split(nfile,"/")[end]
city = split(city,"-")[1]

mkpath(logdir)
mkpath(resdir)


tinit = Dates.now()
tinits = Dates.format(tinit, "ddmmyy-HHhMM-SS")

io = open("$logdir/run$runid-$tinits.log", "w+")
logger = SimpleLogger(io, Logging.Info)
global_logger(logger)
with_logger(logger) do
    @info("--------------------------
    $(Dates.format(tinit, "yy-mm-dd H:M"))
    runid = $runid
    L = [$l1:$dl:$l2], nsamples = $ns
    --------------------------")
    flush(io)

    g, coords, distmx, weightmx, d = OPC.buildCityNetwork(efile, nfile)
    rcellList = OPC.cellList(coords; cellWidth=100.0)

    for ℓ in collect(l1:dl:l2)
        nremoved = []
        dist = []
        origin = []
        destination = []
        seed = Dates.value(DateTime(Dates.now()))
        OD = OPC.odMatrix(ℓ, rcellList; ns=ns, seed=seed, δ = 0.001)
        for sample=1:ns
            # seed = Dates.value(DateTime(Dates.now()))
            # OD = OPC.odMatrix(ℓ, rcellList; ns=1, seed=seed, δ = 0.001)

            (orig, dest)  = OD[sample]
            push!(origin, orig)
            push!(destination, dest)
            p1 = LLA(coords[orig][1], coords[orig][2], 0.0)
            p2 = LLA(coords[dest][1], coords[dest][2], 0.0)
            push!(dist, distance(p1,p2))
            nrem, gr, rmmx = OPC.crackOptimalPaths(g, orig, dest, weightmx)
            if mod(sample, 100)==0
                @info("ℓ = $ℓ, sample = $sample, $(mean(nremoved))")
                flush(io)
            end
            push!(nremoved, nrem)
        end
        fn = "$resdir/$city-nr-$runid-l-$(Int(round(ℓ))).csv"
        df = DataFrame(orig=origin, dest=destination, ell=dist, nr=nremoved)
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
