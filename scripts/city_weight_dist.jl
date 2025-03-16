using CSV, DataFrames
using Random, Distributions
using HypothesisTests
using Plots, Plots.PlotMeasures
using LaTeXStrings
pyplot()
include("utils.jl")
#----

files = ["./data/boston-edges-4h.csv", "./data/manthattan-edges-4h.csv"]
cities = ["Boston", "Manhattan"]
figlab = ["(a)", "(b)"]
colors = [:1, :1]
# files = ["data/boston-edges-4h.csv"]

ktest = 0
pcity = Dict()
for (f,city, color,fl) in zip(files, cities, colors,figlab)
    p1 = plot(
        size = (450, 600),
        legendfontsize = 9,
        top_margin = 3mm,
        bottom_margin = 3mm,
        left_margin=3mm,
        right_margin=3mm,
        legend= :best,
        fg_legend = :white,
        bg_legend = :transparent,
        legendfont=6,
        dpi = 150,
        framestyle = :box,
        xaxis =([0.017, 0.6]),
        yaxis = ([0.21, 5.835]),
        # yaxis = ("\$ \\mathcal{P} \$",font(14), [0., 1]),
        grid = false)
        xaxis!(L"\tau",font(14))
        yaxis!("\$ \\mathcal{P} \$",font(14))

    df = CSV.read(f)|> DataFrame
    n = size(df,1)
    println("n = $n")
    df.tau = df.tt ./ df.len
    n = size(df.tau,1)
    δx = 0.05
    edges = collect(0:δx:0.6)
    pdfd, x = hist(df.tau, edges,norm=true)
    cdfd = δx*cumsum(pdfd)

    dist = fit_mle(LogNormal, df.tau)
    println(dist)
    δxf = 0.01
    xf = 0:δxf:0.6 |> collect
    # xf = x
    pdff = pdf.(dist, xf)
    mf  = round(mean(dist), digits=2)
    sf  = round(std(dist), digits=2)
    flush(stdout); println("$city $mf, $sf")
    cdff = δxf*cumsum(pdff)
    # cdff = cdf.(dist, x)

    yf = Random.rand(dist,n)
    global ktest = ApproximateOneSampleKSTest(df.tau, dist)

    xx = pdfd
    yy = pdf.(dist, x)

    # xx = cdfd
    # yy = δx*cumsum( pdf.(dist, x) )

    mx = mean(xx)
    my = mean(yy)
    vx = sum((xx .- mx).^2)
    vy = sum((yy .- my).^2)
    nn = length(xx)
    Exy = sum( (xx .- mx) .* (yy .- my)  )
    R = Exy/(sqrt(vx)*sqrt(vy))
    R2  = round(R^2, digits=3)
    println("$f R2 = $(R2)")

    p1=bar!(x, pdfd, bar_width=δx, color=color, lab=false, alpha = 0.5)
    p1=plot!(xf, pdff,color=:black, lw=1,lab=false)
    annotate!(0.59, 5.5, text(city, font(12), :right))
    annotate!(0.02, 5.5, text(fl, font(12), :left))

    pcity[city] = p1
    Dn = sqrt(ktest.n)*ktest.δ
    println(Dn)
    B = Kolmogorov()
    Pr = cdf.(B, Dn)
    println(Pr)
end

p2 = plot(pcity["Boston"], pcity["Manhattan"], titlelocation = :left, layout=(2,1))
# savefig(p2, "fig3.pdf")
