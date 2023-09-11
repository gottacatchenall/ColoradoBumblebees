using DrWatson
@quickactivate :ColoradoBumblebees
using CairoMakie

f = Figure(resolution=(3000,2000))

models = [LogisticRegression, BoostedRegressionTree, RandomForest, XGBoost, Ensemble]

treatments = filter(x->length(x)>0,collect(powerset(BEST_REPRESENTATIONS)))
treatments = treatments[sortperm([string(t) for t in treatments])]


mr_dir = joinpath(artifactdir(), "classification_fits", "multiple_representations")
model_dirs = readdir(mr_dir)

rep_dirs = readdir(joinpath(mr_dir, model_dirs[1]))

idx_rep_sorted = sortperm(length.([[x[1] for x in split(r, "_")] for r in rep_dirs]))

rep_dirs = rep_dirs[idx_rep_sorted]


fits = [[ColoradoBumblebees.load(joinpath(mr_dir,m,r)) for m in model_dirs] for r in rep_dirs]


ColoradoBumblebees.load(joinpath(mr_dir, model_dirs[1], rep_dirs[1]))