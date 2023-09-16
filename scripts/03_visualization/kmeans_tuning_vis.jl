using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes
CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(fontsize=42)
set_theme!(fontsize_theme)

dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "KMeans"), readdir(joinpath(artifactdir(), "classification_fits")))]


reps = ColoradoBumblebees.load.(dirs)


K_values = [r[1].embed_model.k for r in representation.(reps)]
prs = praucs.(reps)
rocs = rocaucs.(reps)

category_labels = vcat([[k for _ in eachindex(prs[1])] for k in K_values]...)
markeralpha = 0.5


cols = vcat([[(get(ColorSchemes.colorschemes[:viridis], k/length(K_values)), markeralpha) for _ in eachindex(prs[1])] for k in K_values]...)


I = sortperm(category_labels)

pr_data_array = vcat(prs...)
roc_data_array = vcat(rocs...)
colors = Makie.wong_colors()

pltsettings = (;clouds=nothing,
side=:right,
boxplot_width=0.2,
markersize= 9,
jitter_width=0.2,
boxplot_nudge = 0.5)

begin 
    f = Figure(resolution=(2500, 1200))
    ax = Axis(
        f[1,1],
        xlabel = "k",
        xticks=1:10,
        ylabel = "ROCAUC",
        yticks=0.4:0.1:1
    )
    limits!(ax, 0, 11, 0.4, 1)
    rainclouds!(
        ax, 
        category_labels[I], 
        roc_data_array[I];
        color=cols[I],
        pltsettings...
    )
    hlines!(ax, [0.5], linestyle=:dash, color=:grey50, linewidth=5)

    ax2 = Axis(
        f[1,2],
        xlabel = "k",
        xticks=1:10,
        ylabel = "PRAUC",
        yticks=0.4:0.1:1,)
    limits!(ax2, 0, 11, 0.4, 1)
    rainclouds!(ax2, 
        category_labels[I], 
        pr_data_array[I];    
        color=cols[I], 
    pltsettings...)
    hlines!(ax2, [0.5], linestyle=:dash, color=:grey50, linewidth=5)

    f
end 

save(plotsdir("S1_kmeans.png"), f)