# List of figures:

# Main text:[]
# ----------------------------------------------------------------------
# - F007_predicted_metaweb.png

using DrWatson
using CairoMakie
@quickactivate :ColoradoBumblebees

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=18)
set_theme!(fontsize_theme)


data = load_data()
bee_species, plants_species = bees(data), plants(data)
bee_species, plants_species = bee_species[sortperm([b.name for b in bee_species])], plants_species[sortperm([p.name for p in plants_species])]

binary_prediction, probability_prediction, empirical = ColoradoBumblebees.get_metaweb(BEST_FIT_DIR)

#dir = "/home/michael/Papers/ColoradoBumblebees/artifacts/classification_fits/multiple_representations/ensemble/Temporal"
#binary_prediction, probability_prediction, empirical = ColoradoBumblebees.get_metaweb(dir)

meta = zeros(length(bee_species), length(plants_species))

for (i,b) in enumerate(bee_species), (j,p) in enumerate(plants_species)
    predicted = binary_prediction[b,p]
    observed = empirical[b,p]

    if observed
        meta[i,j] = 1
    end
    if predicted
        #meta[i,j] = 2
    end
    if predicted && observed
        meta[i,j] = 3
    end
end

M = BipartiteNetwork( 
    Matrix{Bool}(meta .> 0), 
    [string(b.name) for b in bee_species], 
    [string(p.name) for p in plants_species]
)

# Make nested 
row_perm = sortperm(sum(eachrow(empirical.adjacency_matrix)))
col_perm = sortperm(sum(eachcol(empirical.adjacency_matrix)))
sorted_bees = empirical.bees[col_perm]
sorted_plants = empirical.plants[row_perm]

meta = meta[col_perm,row_perm]

#E = empirical[col_perm, row_perm]
#P = binary_prediction[col_perm, row_perm]


f = Figure()
ax = Axis(f[1,1])
for i in findall(isequal(1), meta)
    scatter!(ax, [i[1]], [i[2]], color=:green) # only obs
end
for i in findall(isequal(2), meta)
    scatter!(ax, [i[1]], [i[2]], color=:blue) # only predicted
end
for i in findall(isequal(3), meta)
    scatter!(ax,  [i[1]], [i[2]], color=:red)
end
f


bee_label = M.T[col_perm]
plant_label = M.B[row_perm]

mat = meta

begin 

f = Figure(; resolution=(2900, 550))
ax = Axis(
    f[1, 1];
    xticks=(1:(length(plant_label) .+ 0.5), plant_label),
    xticklabelrotation=π / 2,
    yticks=(1:length(bee_label), bee_label),
    xminorgridvisible=true,
    xminorgridcolor=:darkgrey,
    yminorgridcolor=:darkgrey,
    xminorticks=1:(length(plant_label) .+ 0.5),
    yminorticks=1:(length(bee_label) .+ 0.5),
)
limits!(ax, 0.5, length(plant_label) + 0.5, 0.5, length(bee_label) + 0.5)

getcol(val) = [:lightgrey, :red, :deepskyblue2, :teal, ][Int32(val + 1)]

hm = heatmap!(ax, zeros(size(mat))'; colormap=[:grey95])
for (i, b) in enumerate(bee_label), (j, p) in enumerate(plant_label)
    mat[i, j] > 0 &&
        scatter!(ax, [size(mat,2) - j], [i]; marker=:rect, color=getcol(mat[i, j]), markersize=20)
end
f

end

save(plotsdir("F007_predicted_metaweb.png"), f)
save(plotsdir("F007_predicted_metaweb.svg"), f)


# Pt 2:
# Predicted interaction table sorted by probability

df = DataFrame(bee=[], plant=[], probability=[])

predicted_idx = findall(!iszero, P)
for I in predicted_idx
    x,y = I[1], I[2]
    b, p =  sorted_bees[x], sorted_plants[y]
    prob = sorted_probability_prediction[x,y]
    push!(df.bee, b)
    push!(df.plant, p)
    push!(df.probability, prob)
end 

sort!(df, :probability, rev=true)

CSV.write(joinpath(artifactdir(), "predicted_interactions.csv"), df)