using CSV, DataFrames
using Plots, PlotThemes, Plots.PlotMeasures
using Random, Distributions
using LaTeXStrings
pyplot()
theme(:wong2)
PyPlot.rc("text", usetex = "true")
PyPlot.rc("font", family = "CMU Serif")
##

include("utils.jl")
# files = findFilesMaching(r"p0.60-boston-nr-00-l-.+","runs/square/p0.60-boston/")
files = findFilesMaching(r"boston-nr-000-l-.+","runs/cities/boston/")
p1 = plot(
    size = (400, 300),
    legendfontsize = 5,
    legend=:bottomleft,
    top_margin = 3mm,
    bottom_margin = 3mm,
    left_margin=3mm,
    right_margin=3mm,
    fg_legend = :white,
    dpi = 150,
    framestyle = :box,
    xaxis = (L"n_r", font(14)), yaxis=(L"\log_{10}\;P(n_r)", font(14)),
    grid = false,
)
for f in files
    fname = split(f, ".csv")[1]
    fname = split(fname, "/")[end]
    fname = split(fname, "-")[end]

    df = CSV.read(f)|> DataFrame
    nr = df.nr
    n = size(nr,1)
    mi = minimum(nr)
    ma = maximum(nr)

    edges = collect(range(1,stop=500, length=50))
    # edges = collect(range(mi, stop=ma, length=30))
    w = edges[2:end] - edges[1:end-1]
    println("size = $(size(nr,1))")
    h, c = hist(nr, edges, norm=true)
    i = findall(x->x>0.0, h)
    h = h[i]
    c = c[i]
    i = findall(x->x>0.0, c)
    h = h[i]
    c = c[i]
    #l=L"L = "*"$(fname)"*L", \left<nr\right> = "*"$(round.(mean(nr); digits=1))"
    l="\$ L = $(fname) \\;  \\left<nr\\right> = $(round.(mean(nr); digits=1)) \$"
    plot!(c, log10.(h), label=l)
end
annotate!(150, -1.5, Plots.text("Boston road network",7, :right))
plot(p1)
# savefig(p1, "Report20200507/nr_dist_boston.png")
