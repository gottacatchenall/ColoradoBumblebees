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


trait_σ = _splat([representation(r)[1].embed_model.variance_distribution["untruncated"]["σ"] for r in reps], n_reps)
embed_dims = _splat([representation(r)[1].embed_model.truncated_dims for r in reps], n_reps)
numtraits = _splat([representation(r)[1].embed_model.numtraits for r in reps], n_reps)

