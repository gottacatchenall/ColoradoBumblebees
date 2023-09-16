# List of figures:

# Supplement:
# ----------------------------------------------------------------------
# - S005_phylo_node2vec_pr.png   
# - S006_phylo_node2vec_roc.png 


using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes
CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(fontsize=32)
set_theme!(fontsize_theme)

dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "PhylogeneticNode2Vec"), readdir(joinpath(artifactdir(), "classification_fits")))]


reps = ColoradoBumblebees.load.(dirs)

_splat(x, n_reps) = vcat([[i for _ in 1:n_reps] for i in x]...)


prs = vcat(praucs.(reps)...)
rocs = vcat(rocaucs.(reps)...)
n_reps = 128


walklengths = _splat([representation(r)[1].embed_model.walk_length for r in reps], n_reps)
embedding_dim = _splat([representation(r)[1].embed_model.embedding_dim for r in reps], n_reps)
number_of_walks = _splat([representation(r)[1].embed_model.number_of_walks for r in reps], n_reps)



# Its going to be a 4 x 3 panel grid, each with four categories.

# Number of walks: columns (50, 100, 250, 500)
# Walk length: rows (10, 50, 100)
# Embedding dim: elements (4,8,12,16)

axargs = (; titlesize=28,xticks=4:4:16, yticks=0.4:0.05:1)

begin 
fig_roc = Figure(resolution=(2400, 1500))
fig_pr = Figure(resolution=(2400, 1500))

for x in [fig_roc, fig_pr]
    for (i, walklength) in enumerate([10,50,100])
        ax = Axis(x[i,1],
            topspinevisible=false,
            leftspinevisible=false,
            rightspinevisible=false,
            bottomspinevisible=false,
        )
        limits!(ax, 0,1,0,1)
        hidedecorations!(ax)
        text!(ax, 0.2, 0.5, font=:bold, fontsize=36, text="Walk length: $walklength", justification=:right)
    end
end


pltsettings = (;clouds=nothing,
side=:right,
boxplot_width=1,
markersize= 10,
jitter_width=0.5,
boxplot_nudge = 1.5)

markeralpha = 0.3
cols = [(get(ColorSchemes.colorschemes[:viridis], i/6), markeralpha) for i in 1:4]



for (icol, nwalks) in enumerate([50, 100, 250, 500])
    for (irow, walklength) in enumerate([10, 50, 100])
        title = irow == 1 ? "Number of walks: $nwalks" : ""
        xlab = irow == 3 ? "Embedding dimensions" : ""
        ylabpr = icol == 1 ? "PRAUC" : ""
        ylabroc = icol == 1 ? "ROCAUC" : ""
        
        pr_ax = Axis(fig_pr[irow,icol+1]; title=title, ylabel=ylabpr, xlabel=xlab, axargs...)
        limits!(pr_ax, 2.5, 17, 0.4,0.7)

        roc_ax = Axis(fig_roc[irow,icol+1]; title=title, ylabel=ylabroc, xlabel=xlab, axargs...)
        limits!(roc_ax, 2.5,17,0.7, 0.9)

        Iwl = findall(isequal(walklength), walklengths)
        Inw = findall(isequal(nwalks), number_of_walks[Iwl])
        Isorted = sortperm(embedding_dim[Iwl][Inw])

        category_labels = embedding_dim[Iwl][Inw][Isorted]
        
        col = vcat([[cols[findfirst(isequal(i), [4,8,12,16])] for i in category_labels]]...)

        pr_data = prs[Iwl][Inw][Isorted]
        roc_data = rocs[Iwl][Inw][Isorted]

        rainclouds!(pr_ax, category_labels, pr_data; color=col, pltsettings...)
        rainclouds!(roc_ax, category_labels, roc_data; color=col, pltsettings...)

    end
end

fig_pr
fig_roc
end


save(plotsdir("S005_phylo_node2vec_pr.png"), fig_pr)
save(plotsdir("S005_phylo_node2vec_pr.svg"), fig_pr)

save(plotsdir("S006_phylo_node2vec_roc.png"), fig_roc)
save(plotsdir("S006_phylo_node2vec_roc.svg"), fig_roc)