using DataFrames, CSV

df = CSV.read(datadir("public", "occurrence", "bees.csv"), DataFrame)


onespecies = filter(x->x.species==unique(df[!,:species])[1], df)

long, lat = onespecies[!, :longitude],onespecies[!,:latitude]

using CairoMakie
using Clustering


X = hcat(long, lat)'
K = 5
res = kmeans(X, K)

ctrs = [res.centers[:,i] for i in 1:K]


fig = Figure()
ax1 = Axis(fig[1,1])
scatter!(ax1, long,lat, color=(:blue, 0.05))
scatter!(ax1, [c[1] for c in ctrs], [c[2] for c in ctrs], markersize=20, color=:orange)
current_figure()