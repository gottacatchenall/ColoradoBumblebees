# List of figures:

# Supplement:
# ----------------------------------------------------------------------
# - S017_metaweb_svd_tuning_pr.png   
# - S018_metaweb_svd_tuning_roc.png 


using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes
CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(fontsize=32)
set_theme!(fontsize_theme)
n_reps =128

svd_dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "MetawebSVD"), readdir(joinpath(artifactdir(), "classification_fits")))]
lfsvd_dirs = [joinpath(artifactdir(), "classification_fits", x) for x in filter(x->contains(x, "LFSVD"), readdir(joinpath(artifactdir(), "classification_fits")))]

svd_reps = ColoradoBumblebees.load.(svd_dirs)
lfsvd_reps = ColoradoBumblebees.load.(lfsvd_dirs)

prs_svd = vcat(praucs.(svd_reps)...)
prs_lfsvd = vcat(praucs.(lfsvd_reps)...)


rocs_svd =  vcat(rocaucs.(svd_reps)...)
rocs_lfsvd = vcat(rocaucs.(lfsvd_reps)...)

trunc_dims = sort(unique([representation(r)[1].embed_model.truncation_dims for r in svd_reps]))
embed_dims = sort(unique([representation(r)[1].embed_model.embed_dims for r in svd_reps]))

truncs = vcat([[representation(r)[1].embed_model.truncation_dims for _ in 1:n_reps] for r in svd_reps]...)
embed = vcat([[representation(r)[1].embed_model.embed_dims for _ in 1:n_reps] for r in svd_reps]...)



category_labels = []
svd_pr_data_matrix = []
svd_roc_data_matrix = []
lfsvd_pr_data_matrix = []
lfsvd_roc_data_matrix = []



for e in embed_dims
    for t in trunc_dims
        It = findall(isequal(t), truncs)
        Ie = findall(isequal(e), embed[It])
        

        push!(category_labels, [("($e, $t)") for _ in 1:length(prs_svd[It][Ie])])
        push!(svd_pr_data_matrix, prs_svd[It][Ie])
        push!(svd_roc_data_matrix, rocs_svd[It][Ie])
        push!(lfsvd_pr_data_matrix, prs_lfsvd[It][Ie])
        push!(lfsvd_roc_data_matrix, rocs_lfsvd[It][Ie])

    end
end

category_labels = vcat(category_labels...)
svd_pr_data_matrix = vcat(svd_pr_data_matrix...)
svd_roc_data_matrix = vcat(svd_roc_data_matrix...)
lfsvd_pr_data_matrix = vcat(lfsvd_pr_data_matrix...)
lfsvd_roc_data_matrix = vcat(lfsvd_roc_data_matrix...)


pltsettings = (;clouds=nothing,
    side=:right,
    boxplot_width=0.2,
    markersize= 10,
    jitter_width=0.2,
    boxplot_nudge = -0.1)

markeralpha = 0.5

cat_idx = map(x->findfirst(isequal(x), unique(category_labels)), category_labels)

cols = [(get(ColorSchemes.colorschemes[:viridis], i/20), markeralpha) for i in cat_idx]
    

begin
    fig_pr = Figure(resolution=(2400, 1400))

    ax_svd = Axis(fig_pr[1,1], title="SVD", xticks=0.5:0.1:1)
    ylims!(ax_svd, 0.7, 0.91)
    ax_lfsvd = Axis(fig_pr[2,1], title="Linear-Filter SVD", xlabel = "(Embedding Dims, Truncate Dims)")
    ylims!(ax_lfsvd, 0.7, 0.91)

    rainclouds!(ax_svd, category_labels, svd_pr_data_matrix; color=cols, ylabel="PRAUC", pltsettings...)
    rainclouds!(ax_lfsvd, category_labels, lfsvd_pr_data_matrix; color=cols, ylabel="PRAUC", pltsettings...)

    fig_pr
end


begin 
    fig_roc = Figure(resolution=(2400, 1400))
    ax_svd = Axis(fig_roc[1,1], title="SVD", yticks=0.8:0.025:1, ylabel="ROCAUC")
    ax_lfsvd = Axis(fig_roc[2,1], title="Linear-Filter SVD",yticks=0.8:0.025:1, xlabel = "(Embedding Dims, Truncate Dims)", ylabel="ROCAUC")
    ylims!(ax_svd, 0.875, 0.95)
    ylims!(ax_lfsvd, 0.875, 0.95)

    rainclouds!(ax_svd, category_labels, svd_roc_data_matrix; color=cols, pltsettings...)
    rainclouds!(ax_lfsvd, category_labels, lfsvd_roc_data_matrix; color=cols, pltsettings...)

    fig_roc
end 

save(plotsdir("S017_metaweb_svd_tuning_pr.png"), fig_pr)
save(plotsdir("S017_metaweb_svd_tuning_pr.svg"), fig_pr)

save(plotsdir("S018_metaweb_svd_tuning_roc.png"), fig_roc )
save(plotsdir("S018_metaweb_svd_tuning_roc.svg"), fig_roc )