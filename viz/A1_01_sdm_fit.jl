"""
    This script plots the ROC-AUC and TSS for each SDM, split by bees and flowers.
"""

using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG

include(joinpath("..", "src", "io.jl"))

CairoMakie.activate!(; px_per_unit=3)
sdms = read_sdms("./artifacts")


Dict([k=>v[:metrics]["mcc"] for (k,v) in sdms]) |> findmin

species = collect(keys(sdms))

plant_color = (colorant"#5e9e9dff", 0.2)
bee_color = (colorant"#5f95c7ff", 0.4)



X = [occursin("Bombus", x) ? 1 : 1.2 for x in species]
cols = [x == 1 ? bee_color : plant_color for x in X]
μ_mcc = [sdms[sp][:metrics]["mcc"] for sp in species]



# Plot

begin 
f = Figure()
shared_axis_args = (;
    xticks=([1, 1.2], ["Bumble bees", "Plants"]), 
    xgridvisible = false,
    xticksvisible = false,
    xticklabelsize = 25, 
    ylabelsize= 25,
    yticklabelsize = 18
)
shared_raincloud_args = (;
    color=cols, 
    gap=0.03, 
    side_nudge = 0.02, 
    markersize=17, 
    cloud_width=0.1, 
    boxplot_width=0.03, 
    jitter_width=0.02
)
ax2 = Axis(
    f[1,1];
    ylabel="MCC",
    shared_axis_args...
)
limits!(ax2, 0.85, 1.3, 0.4, 1.)
rc2 = rainclouds!(ax2, X, μ_mcc; shared_raincloud_args...)
f
end 

mean(μ_mcc)

isbee = [occursin("Bombus", x) for x in species]

μ_mcc[isbee] |> median
μ_mcc[.!isbee] |> median


save("plots/sdm_fit.png", f)

