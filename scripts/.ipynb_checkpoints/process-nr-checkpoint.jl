using CSV, DataFrames
using Plots, Plots.PlotMeasures
using Statistics
using PlotThemes
using LaTeXStrings
using LsqFit
pyplot()
theme(:wong)
PyPlot.rc("text", usetex = "true")
PyPlot.rc("font", family = "CMU Serif")
include("utils.jl")
##
mks = [:circle,:utriangle, :cross,  :star5,
        :diamond,  :rect, :hexagon,  :xcross,
        :dtriangle, :rtriangle, :ltriangle,
        :pentagon, :heptagon, :octagon,
        :star4, :star6, :star7, :star8, :vline, :hline,  :x, :+]


joinfiles = false

# folders=["runs/square/p0.60-beta-0.002",
#          "runs/square/p0.60-beta-1.800",
#          "runs/square/p0.60-beta-200.0",
#          "runs/square/p0.60-boston"]
# labels=[L"\beta = 0.002",
#         L"\beta = 1.800",
#         L"\beta = 200.0",
#         "Boston\nweight distrib."]


# folders=["runs/square/p0.00-boston/",
#          "runs/square/p0.10-boston/",
#          "runs/square/p0.20-boston/",
#          "runs/square/p0.30-boston/",
#          "runs/square/p0.40-boston/",
#          "runs/square/p0.50-boston",
#          "runs/square/p0.60-boston",
#          "runs/square/p0.70-boston",
#          "runs/square/p0.80-boston",
#          "runs/square/p0.90-boston",
#          "runs/square/p0.95-boston"]
# labels=[L"p = 0.0",
#         L"p = 0.1",
#         L"p = 0.2",
#         L"p = 0.3",
#         L"p = 0.4",
#         L"p = 0.5",
#         L"p = 0.6",
#         L"p = 0.7",
#         L"p = 0.8",
#         L"p = 0.9",
#         L"p = 0.95"]

# compara square network com boston network
folders=["runs/square/p0.70-boston/",
         "runs/square/p0.74-boston",
         "runs/square/p0.80-boston",
         "runs/cities/boston"
         ]
labels=["square "*L"p = 0.70",
        "square "*L"p = 0.74",
        "square "*L"p = 0.80",
        "Boston network\n(p=0.4)"]


p1 = plot(
    size = (400, 300),
    legendfontsize = 9,
    top_margin = 3mm,
    bottom_margin = 3mm,
    left_margin=3mm,
    right_margin=3mm,
    legend= bbox_to_anchor = (0.72, 0.45),
    fg_legend = :white,
    bg_legend = :white,
    legendfont=6,
    # legendtitle= L"p = 0.6",
    legendtitlefont= 7,
    dpi = 150,
    framestyle = :box,
    xaxis =(L"\log_{10}(L)",font(14)),
    yaxis = (L"\log_{10}(N_r)",font(14)),
    grid = false)

@. model(x, p) = p[1]*x + p[2]

for (i,(folder, lab)) in enumerate(zip(folders,labels))
    city = split(folder,"/")[end]
    println(city)
    lav = []
    lerr = []
    nrav = []
    nrerr = []
    files = findFilesMaching("$(city)"r"-nr-.+-l-.+\.csv", folder)
    L = []
    for f in files
        f = split(f,"-l-")[end]
        f = split(f,r"[.,-]")[1]
        push!(L, f)
    end
    L = map(x->parse(Int,x),L)
    sort!(L)
    unique!(L)
    for (j,ℓ) in enumerate(L)
        df = DataFrame()
        files = findFilesMaching("$(city)"r"-nr-.+-l-"*"$(ℓ)"r".csv", folder)
        for f in files
            dfi = CSV.read(f)
            df = [df; dfi]
        end
        df = unique(df)
        if joinfiles
            fn = "$folder/c-$city-nr-000-l-$(ℓ).csv"
            CSV.write(fn, df)
        end
        n = size(df,1)
        if j==1
            println("$folder $(size(files,1)) files, $n points")
        end
        d = mean(df.ell)
        derr = std(df.ell)/sqrt(n)
        nr = mean(df.nr)
        nre = std(df.nr)/sqrt(n)
        push!(lav,d)
        push!(lerr, derr)
        push!(nrav, nr)
        push!(nrerr, nre)
    end
    x = log10.(lav)
    y = log10.(nrav)

    # x = lav
    # y = nrav
    ms = 8
    if i==4
        ms=10
    end
    p1 = scatter!(x,y,
                 markershape=mks[i], markersize=ms,
                 label=lab, color=i)
    p0 = [0.5, 0.5]
    fit = curve_fit(model, x,y,p0)
    yf = model(x, fit.param)
    println(fit.param)
    plot!(x,yf, lw=0.5, ls=:dash, color=i, label="$( round(fit.param[1], digits=2) )")
end
# p1=annotate!(3.25, 2.9, Plots.text("all  using  Boston weight distribution",7, :left))
# p1=annotate!(3.25, 2.47, Plots.text("p = 0.6",7, :left))

savefig(p1,"Report20200507/Nr_boston.png")
# savefig(p1,"Report20200507/p06.png")
p1
