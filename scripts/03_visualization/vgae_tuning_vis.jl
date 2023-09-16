using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes
CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(fontsize=32)
set_theme!(fontsize_theme)

dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "GraphAutoencoder"), readdir(joinpath(artifactdir(), "classification_fits")))]

variational_dirs = filter(x->contains(x, "Variational"), dirs)
standard_dirs = filter(x->contains(x, "Standard"), dirs)

n_reps = 128

standard_reps = ColoradoBumblebees.load.(standard_dirs)
variational_reps = ColoradoBumblebees.load.(variational_dirs)

standard_pr = vcat(praucs.(standard_reps)...)
variational_pr = vcat(praucs.(variational_reps)...)

standard_roc = vcat(rocaucs.(standard_reps)...)
variational_roc = vcat(rocaucs.(variational_reps)...)


input_features = sort(unique([representation(r)[1].embed_model.num_input_features for r in standard_reps]))
embed_dims = sort(unique([representation(r)[1].embed_model.embedding_dim for r in standard_reps]))

_splat(x, n_reps) = vcat([[i for _ in 1:n_reps] for i in x]...)

num_input_features = _splat([representation(r)[1].embed_model.num_input_features for r in standard_reps], n_reps)
embedding_dim = _splat([representation(r)[1].embed_model.embedding_dim for r in standard_reps], n_reps)

# Each col is for each input feature num (16,64,256,1028)
# Each panel has 5 elements per embedding dims (4,8,16,32,64)
# Top row is Standard, Bottom row is Variational

pltsettings = (;clouds=nothing,
side=:right,
boxplot_width=0.25,
markersize= 10,
jitter_width=0.25,
boxplot_nudge = 0.6)

cols = [(get(ColorSchemes.colorschemes[:viridis], i/8), markeralpha) for i in 1:5]


begin 

    fig_pr = Figure(resolution=(2400,1400))
    
    ax1 = Axis(fig_pr[1,1],
        topspinevisible=false,
        leftspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
    )
    ax2 = Axis(fig_pr[2,1],
    topspinevisible=false,
    leftspinevisible=false,
    rightspinevisible=false,
    bottomspinevisible=false)
    hidedecorations!(ax1)
    hidedecorations!(ax2)
    xlims!(ax1, 0,1)
    xlims!(ax2, 0,1)
    text!(ax1, 0.2, 0.5, font=:bold, fontsize=36, text="Standard GAE", justification=:right)
    text!(ax2, 0.2, 0.5, font=:bold, fontsize=36, text="Variational GAE", justification=:right)

    for (i, input_feat) in enumerate(input_features)
        ax = Axis(fig_pr[1,i+1],
            xticksvisible=false,
            xticklabelsvisible=false,
            ylabel= i == 1 ? "PRAUC" : ""
        )
        ylims!(ax, 0.4, 0.75)
        ax2 = Axis(
            fig_pr[2,i+1],
            xlabel="Embedding Dimensions",
            ylabel = i == 1 ? "PRAUC" : "",
            xticks = (1:5, string.(embed_dims))
        )
        ylims!(ax2, 0.35, 0.75)

        Inif = findall(isequal(input_feat), num_input_features)
        category_labs = map(x->findfirst(isequal(x), embed_dims), embedding_dim[Inif])
        
        col = [cols[i] for i in category_labs]
        
        rainclouds!(ax, category_labs, standard_pr[Inif]; 
            title = "Input Features: $input_feat",        
            color=col, pltsettings...)
        rainclouds!(ax2, category_labs, variational_pr[Inif]; color=col, pltsettings...)

    end

    fig_pr

end 

begin 

    fig_roc = Figure(resolution=(2400,1400))
    ax1 = Axis(fig_roc[1,1],
    topspinevisible=false,
    leftspinevisible=false,
    rightspinevisible=false,
    bottomspinevisible=false,
    )
    ax2 = Axis(fig_roc[2,1],
        xlabel="Embedding Dimensions",
        topspinevisible=false,
        leftspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
        xticks = (1:5, string.(embed_dims))
    )
    hidedecorations!(ax1)
    hidedecorations!(ax2)
    xlims!(ax1, 0,1)
    xlims!(ax2, 0,1)
    text!(ax1, 0.2, 0.5, font=:bold, fontsize=36, text="Standard GAE", justification=:right)
    text!(ax2, 0.2, 0.5, font=:bold, fontsize=36, text="Variational GAE", justification=:right)

    for (i, input_feat) in enumerate(input_features)
        ylab = 
        ax = Axis(fig_roc[1,i+1], 
            ylabel=i == 1 ? "ROCAUC" : "",
            xticksvisible=false,
            xticklabelsvisible=false
        )
        ylims!(ax, 0.725, 0.875)
        ax2 = Axis(fig_roc[2,i+1], 
        ylabel = i == 1 ? "ROCAUC" : "",            
        xticks = (1:5, string.(embed_dims)),
        xlabel="Embedding Dimensions",

        )
        ylims!(ax2, 0.725, 0.875)

        Inif = findall(isequal(input_feat), num_input_features)
        category_labs = map(x->findfirst(isequal(x), embed_dims), embedding_dim[Inif])
        col = [cols[i] for i in category_labs]
        rainclouds!(ax, category_labs, standard_roc[Inif];             
            title = "Input Features: $input_feat",   
            color=col, pltsettings...)
        rainclouds!(ax2, category_labs, variational_roc[Inif]; color=col, pltsettings...)

    end

    fig_roc
end 


save(plotsdir("S7_metaweb_gae_pr.png"), fig_pr)
save(plotsdir("S8_metaweb_gae_roc.png"), fig_roc)