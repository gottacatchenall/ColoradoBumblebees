# List of figures:

# Main text:
# ----------------------------------------------------------------------
# - F003_representation_comparison.png
#
# Supplement
# - S019_ROC_PR_correlation.png

using DrWatson
@quickactivate :ColoradoBumblebees
using CairoMakie


CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=32)
set_theme!(fontsize_theme)


f = Figure(resolution=(3000,2000))

models = [LogisticRegression, BoostedRegressionTree, RandomForest, XGBoost, Ensemble]

treatments = filter(x->length(x)>0,collect(powerset(BEST_REPRESENTATIONS)))
treatments = treatments[sortperm([string(t) for t in treatments])]


mr_dir = joinpath(artifactdir(), "classification_fits", "multiple_representations")
model_dirs = readdir(mr_dir)

rep_dirs = readdir(joinpath(mr_dir, model_dirs[1]))

model_dirs = model_dirs[[3,4,1,5,2]]

idx_rep_sorted = sortperm(length.([[x[1] for x in split(r, "_")] for r in rep_dirs]))
rep_dirs = rep_dirs[idx_rep_sorted]


fits = [[ColoradoBumblebees.load(joinpath(mr_dir,m,r)) for m in model_dirs] for r in rep_dirs]

model_ord = ["Logistic", "Random Forest", "Boosted Regression Tree", "XGBoost", "Ensemble"]



prs = [praucs.(f) for f in fits]
rocs = [rocaucs.(f) for f in fits]

roc_mat = zeros(length(rep_dirs), length(model_dirs))
pr_mat = zeros(length(rep_dirs), length(model_dirs))


for (i,r) in enumerate(rep_dirs)
    for (j,model) in enumerate(model_dirs)
        pr_mat[i,j] = mean(prs[i][j])
        roc_mat[i,j] = mean(rocs[i][j])
    end
end

hmsettings = (; colorrange=(0.5,1), colormap=:thermal)

begin

f = Figure(resolution=(2000,1000))

ax = Axis(
    f[1,1], 
    title="Precision-Recall Area-Under-the-Curve (PR-AUC)",
    yticks=(1:5, model_ord),
    yticksvisible=false,
    titlealign=:left
)

hm1 = heatmap!(ax, pr_mat; hmsettings...)

for (i,r) in enumerate(rep_dirs)
    for (j,model) in enumerate(model_dirs)
        text!(ax, i-0.3,j-0.25, text="$(string(mean(prs[i][j]))[2:4])", color=:white, fontsize=20)
    end
end

ax2 = Axis(f[2,1];
    title="Receiver-Operator-Curve Area-Under-the-Curve (ROC-AUC)",
    titlealign=:left,
    yticks=(1:5, model_ord),
    yticksvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false
)

hm2 = heatmap!(ax2, roc_mat; hmsettings...)

for (i,r) in enumerate(rep_dirs)
    for (j,model) in enumerate(model_dirs)
        text!(ax2, i-0.3,j-0.25, text="$(string(mean(rocs[i][j]))[2:4])", color=:white, fontsize=20)
    end
end

Colorbar(f[:,2], hm2, width=25, label="AUC")

xlabels = [sort([y[1] for y in split(x, "_")]) for x in rep_dirs]


# Structural: blue,
# Temporal: Navy 
# Relative Abundnace: Purple
# Environment: Orange
# Phylogenetic: Red 
colorvec = Dict("S"=>colorant"#81a8c1", "T"=>colorant"#4c566a", "E"=>colorant"#d08770", "R"=>"#b48EAD", "P"=>colorant"#bf616a")

ticklabs = ["Temporal", "Structural", "Relative Abundance", "Phylogenetic", "Environment"]

xlabel_ax = Axis(
    f[3,:],
    yticksvisible=false,
    yticks=(0:4, ticklabs),
    xticksvisible=false,
    xticklabelsvisible=false
)
xlims!(xlabel_ax, 0,32)
ylims!(xlabel_ax, -0.5,5)
#hidedecorations!(xlabel_ax)
hidespines!(xlabel_ax)

yheight = Dict([y=>5-findfirst(isequal(y), xlabels[end]) for y in xlabels[end]])

for x in 1:length(rep_dirs)
    for (y,char) in enumerate(xlabels[x])
        scatter!(xlabel_ax, x-0.5, yheight[char], color=colorvec[string(char)], marker=:rect, markersize=55)
    end
end

f

save(plotsdir("F003_representation_comparison.png"), f)
save(plotsdir("F003_representation_comparison.svg"), f)

end 





begin
xs = vcat([[mean(i) for i in x] for x in rocs]...)
ys = vcat([[mean(i) for i in x] for x in prs]...)

f = Figure(resolution=(1200,1200))
ax = Axis(f[1,1], xlabel="ROC-AUC", ylabel="PR-AUC", xticks=0.4:0.1:1, yticks=0.4:0.1:1)
limits!(0.4, 1, 0.4 ,1)
scatter!(ax, xs,ys, markersize=30, color=(:dodgerblue, 0.4))
lines!(ax, [0, 1], [0, 1])
f
save(plotsdir("S019_ROC_PR_correlation.png"), f)
save(plotsdir("S019_ROC_PR_correlation.svg"), f)

end