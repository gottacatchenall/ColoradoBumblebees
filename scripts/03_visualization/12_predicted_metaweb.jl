# List of figures:

# Main text:
# ----------------------------------------------------------------------
# - F007_predicted_metaweb.png

using DrWatson
using CairoMakie
@quickactivate :ColoradoBumblebees

data = load_data()
bee_species, plants_species = bees(data), plants(data)
bee_species, plants_species = bee_species[sortperm([b.name for b in bee_species])], plants_species[sortperm([p.name for p in plants_species])]

binary_prediction, probability_prediction, empirical = ColoradoBumblebees.get_metaweb(BEST_FIT_DIR)

M = BipartiteNetwork( Matrix{Bool}(any.(binary_prediction .∪ empirical)), [string(b.name) for b in bee_species], [string(p.name) for p in plants_species],)

# Make nested 
row_perm = sortperm(sum(eachrow(adjacency(M))))
col_perm = sortperm(sum(eachcol(adjacency(M))))

sorted_bees = M.T[col_perm]
sorted_plants = M.B[row_perm]

E = empirical[col_perm, row_perm]
P = binary_prediction[col_perm, row_perm]

f = Figure()
ax = Axis(f[1,1])
for i in findall(E)
    scatter!(ax, [i[1]], [i[2]], color=:blue)
end
for i in findall(P)
    scatter!(ax,  [i[1]], [i[2]], color=:red)
end
f


bee_label = M.T[col_perm]
plant_label = M.B[row_perm]

mat = zeros(Int64,size(E))
mat[findall(!iszero, P)] .= 2
mat[findall(!iszero, E)] .= 1

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

getcol(val) = [:lightgrey, :deepskyblue2, :teal, :red][val + 1]

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

predicted_idx = findall(!iszero, P)
for I in predicted_idx
    @info I
    @info sorted_bees[I[1]], sorted_plants[I[2]]
end 