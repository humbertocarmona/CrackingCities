using OptimumPathCrack
using Geodesy
using Dates
using Statistics
using CSV, DataFrames
using Memento
using ArgParse
using SparseArrays

logger = getlogger("OptimumPathCrack")
# %%  --------------------------------------------------------------------------
# googlekey = "AIzaSyApQzC_OLdxiITS7ynh_XsWZZOU8XOKQHs"

function parse_commandline()
	s = ArgParseSettings();

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
parsed_args = Dict{String, Any}(
	"l1" => 1000.0,
	"l2" => 1000.0,
	"dl" => 1.0,
	"nsamples" => 10,
	"runid" => "1",
	"efile" => "../data/boston-edges-4h.csv",
	"nfile" => "../data/boston-nodes.csv",
	"logdir" => "logs",
	"resdir" => "results",
);

l1 = parsed_args["l1"];
l2 = parsed_args["l2"];
dl = parsed_args["dl"];
n_samples = parsed_args["nsamples"];
runid = lpad(parsed_args["runid"], 3, "0");
efile = parsed_args["efile"];
nfile = parsed_args["nfile"];
logdir = parsed_args["logdir"];
resdir = parsed_args["resdir"];
city = split(nfile, "/")[end];
city = split(city, "-")[1];

info(logger, "creating $(logdir)")
mkpath(logdir);
info(logger, "creating $(resdir)")
mkpath(resdir);

# Create a file handler and attach it
log_file = joinpath(pwd(), logdir, "myLogFile.log");
handler = DefaultHandler(
	log_file,
	DefaultFormatter("{level}: {msg}"),
);
setlevel!(logger, "debug");

push!(logger, handler);

time_begin = Dates.now();
time_begin_str = Dates.format(time_begin, "ddmmyy-HHhMM-SS");
info(logger, "time_begin =  $(time_begin_str)");
info(logger, "runid =  $(runid)");
info(logger, "L =  [$(l1):$(dl):$(l2)]");
info(logger, "nsamples =  $(n_samples)");

# build the graph with travel time and distances
g, coords, distance_matrix, weight_matrix, edges_index_dict = buildCityNetwork(efile, nfile)

# Creates a spatial cell list for "coord" coordinates, cellWidth in meters
cell_list = cellList_lat_lon(coords; cellWidth = 100.0);
# Dict("coords" => coords,
# 		"pos" => pos,
# 		"dx" => d_lon, "dy" => d_lat,
# 		"nx" => nx, "ny" => ny,
# 		"cells" => cells,
# 		"next" => next)


for ℓ in collect(l1:dl:l2)
	info(logger, "ℓ = $(ℓ)");
	nremoved = [];
	dist = [];
	origin = [];
	destination = [];
	seed = Dates.value(DateTime(Dates.now()));
	# generate n_samples OD pairs with distance ℓ ± δ
	OD = create_ODs(ℓ, cell_list, n_od_pairs = n_samples, seed = seed, δ = 0.001);

	for sample ∈ 1:n_samples
		# seed = Dates.value(DateTime(Dates.now()))
		(orig, dest) = OD[sample];

		push!(origin, orig);
		push!(destination, dest);
		p1 = LLA(coords[orig][1], coords[orig][2], 0.0);
		p2 = LLA(coords[dest][1], coords[dest][2], 0.0);
		push!(dist, euclidean_distance(p1, p2));

		gr, removed_edges, path_edges = crackOptimalPaths(g, orig, dest, weight_matrix);
        nrem = length(findnz(removed_edges)[3]);
		push!(nremoved, nrem);

		if mod(sample, 100) == 0
			info(logger, "ℓ = $ℓ, sample = $sample, $(mean(nremoved))");
		end
	end
	fn = "$resdir/$city-nr-$runid-l-$(Int(round(ℓ))).csv"
	df = DataFrame(orig = origin, dest = destination, ell = dist, nr = nremoved);
	CSV.write(fn, df);
	info(logger, "wrote $fn")
end
tend = Dates.now();
dur = Dates.canonicalize(Dates.CompoundPeriod(tend - time_begin));
mess = "----------- finished -----------------
$(Dates.format(tend, "yy-mm-dd H:M"))
took  $dur
--------------------------------------------";
info(logger, mess);

