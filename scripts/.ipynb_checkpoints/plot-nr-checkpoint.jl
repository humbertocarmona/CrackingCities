using Plots, Plots.PlotMeasures
using LaTeXStrings
using LsqFit
using Printf
using CSV
using DataFrames


pyplot()
PyPlot.rc("text", usetex = "true")
PyPlot.rc("font", family = "CMU Serif")

df = CSV.read("results/square-p0.6-beta0.005-500-nr.csv")
dfb = CSV.read("data/boston-nr.csv")
dfm = CSV.read("data/manthattam-nr.csv")
@. model(x, p) = p[2] * x^p[1]

fit = curve_fit(model, collect(df.L), collect(df.nr), [1.0, 1.0])
println(fit.param)
p1 = plot(
    size = (400, 300),
    legendfontsize = 9,
    legend=:topleft,
    top_margin = 3mm,
    bottom_margin = 3mm,
    left_margin=3mm,
    right_margin=3mm,
    fg_legend = :white,
    dpi = 150,
    framestyle = :box,
    xaxis =(L"\ell", :log10, font(15)),
    yaxis = (L"N_r", :log10, font(15)),
    grid = false,
    # thickness_scaling=1.1
)
scatter!(df.L, df.nr,label = L"p = 0.6\;\;\beta = 0.002",
         markersize=5, c=1, yerror=df.nrerr, lw=2)
x=1e3:1e2:1e4; y=model(x, fit.param)
plot!(x, y, label="", c=:black)


scatter!(dfb.L, dfb.nr,label = L"Boston",
         markersize=5, c=2, yerror=dfb.nrerr, lw=2)


scatter!(dfm.L, dfm.nr,label = L"Manhattan",
                  markersize=5, c=3, yerror=dfm.nrerr, lw=2)

# savefig(p1, "Results/nr.png")
