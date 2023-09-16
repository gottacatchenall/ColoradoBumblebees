using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes
CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(fontsize=32)
set_theme!(fontsize_theme)

dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "SimulatedTraits"), readdir(joinpath(artifactdir(), "classification_fits")))]


reps = ColoradoBumblebees.load.(dirs)

prs = vcat(praucs.(reps)...)
rocs = vcat(rocaucs.(reps)...)
n_reps = 128

_splat(x, n_reps) = vcat([[i for _ in 1:n_reps] for i in x]...)


trait_σ = _splat([representation(r)[1].embed_model.stddev_of_stddev for r in reps], n_reps)
embed_dims = _splat([representation(r)[1].embed_model.truncated_dims for r in reps], n_reps)
numtraits = _splat([representation(r)[1].embed_model.numtraits for r in reps], n_reps)


unique_trait_σ = sort(unique(trait_σ))
unique_embed_dims = sort(unique(embed_dims))
unique_num_traits = sort(unique(numtraits))

pltsettings = (;clouds=nothing,
side=:right,
boxplot_width=0.75,
markersize= 9,
jitter_width=0.4,
boxplot_nudge = 1.1)

markeralpha = 0.5
cols = [(get(ColorSchemes.colorschemes[:viridis], i/10), markeralpha) for i in 1:length(unique_embed_dims)]


begin 
f = Figure(resolution=(2000,2000))
ax1 = Axis(f[1,1])
ax2 = Axis(f[2,1])
ax3 = Axis(f[3,1])
hidespines!.([ax1, ax2, ax3])
hidedecorations!.([ax1, ax2, ax3])
xlims!.([ax1, ax2, ax3], -0.1, 0.5)
text!(ax1, 0., 0.5, text="Trait Std. Dev:  $(unique_trait_σ[1])")
text!(ax2, 0., 0.5, text="Trait Std. Dev:  $(unique_trait_σ[2])")
text!(ax3, 0., 0.5, text="Trait Std. Dev: $(unique_trait_σ[3])")

for (i,σ) in enumerate(unique_trait_σ)
    for (j,nt) in enumerate(unique_num_traits)
        xlabel = i == 3 ? "Embed Dims" : ""
        ylabel = j == 1 ? "PRAUC" : ""

        title = i == 1 ? "Number of traits: $nt" : ""
        ax = Axis(
            f[i,j+1];
            subtitle=title,
            xlabel=xlabel,
            ylabel=ylabel,
            xticks=2:2:16,
            yticks=0:0.05:1
        )
        ylims!(ax, 0.4, 0.725)
        xlims!(ax, 0,17)
        I = findall(isequal(σ), trait_σ)
        J = findall(isequal(nt), numtraits[I])

        category_labels = embed_dims[I][J]
        col = vcat([[cols[findfirst(isequal(i), unique_embed_dims)] for i in category_labels]]...)

        y = prs[I][J]
        rainclouds!(ax, category_labels, y; color=col, pltsettings...)
    end
end

f
end
save(plotsdir("S15_phylo_simulated_traits_pr.png"), f)
begin 
    f = Figure(resolution=(2000,2000))
    ax1 = Axis(f[1,1])
    ax2 = Axis(f[2,1])
    ax3 = Axis(f[3,1])
    hidespines!.([ax1, ax2, ax3])
    hidedecorations!.([ax1, ax2, ax3])
    xlims!.([ax1, ax2, ax3], -0.1, 0.5)
    text!(ax1, 0., 0.5, text="Trait Std. Dev:  $(unique_trait_σ[1])")
    text!(ax2, 0., 0.5, text="Trait Std. Dev:  $(unique_trait_σ[2])")
    text!(ax3, 0., 0.5, text="Trait Std. Dev: $(unique_trait_σ[3])")
    
    for (i,σ) in enumerate(unique_trait_σ)
        for (j,nt) in enumerate(unique_num_traits)
            xlabel = i == 3 ? "Embed Dims" : ""
            ylabel = j == 1 ? "ROCAUC" : ""
    
            title = i == 1 ? "Number of traits: $nt" : ""
            ax = Axis(
                f[i,j+1];
                subtitle=title,
                xlabel=xlabel,
                ylabel=ylabel,
                xticks=2:2:16,
                yticks=0:0.05:1
            )
            ylims!(ax, 0.7, 0.9)
            xlims!(ax, 0,17)
            I = findall(isequal(σ), trait_σ)
            J = findall(isequal(nt), numtraits[I])
    
            category_labels = embed_dims[I][J]
            col = vcat([[cols[findfirst(isequal(i), unique_embed_dims)] for i in category_labels]]...)
    
            y = rocs[I][J]
            rainclouds!(ax, category_labels, y; color=col, pltsettings...)
        end
    end
    
    f
    end

    save(plotsdir("S16_phylo_simulated_traits_roc.png"), f)