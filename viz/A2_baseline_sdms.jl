"""
    This script is used to create the baseline SDMs for the appendix.
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
include("shared.jl")


CairoMakie.activate!(; px_per_unit=3)
sdms = read_sdms("./artifacts")

state_poly, county_poly = get_polygons()

function plot_baseline_sdm(
    fig,
    slice,
    sdm_obj;
    title = "",
)
    range = sdm_obj[:baseline][:range]

    ax = Axis(
        fig[slice...],
        titlealign=:left,
        titlefont=:bold_italic,
        aspect = 1.,
        title = title,
    )
    heatmap!(ax, range, colormap=[:grey95, :seagreen4])
    scatter!(ax, Bool.(sdm_obj[:baseline][:presences]), color=:yellow, markersize=2)
    scatter!(ax, Bool.(sdm_obj[:baseline][:absences]), color=:black, markersize=2)
    lines!(ax, county_poly, color=:black, linewidth=0.25)
    lines!(ax, state_poly, color=:black, linewidth=1.)
end


f = Figure()
plot_baseline_sdm(f, (1,1), first(sdms)[2])
f

scatter(Bool.(first(sdms)[2][:baseline][:presences]))





outdir = joinpath("plots", "baseline_sdms")
mkpath(outdir)

spnames = sort(collect(keys(sdms)))
cidx = CartesianIndices((1:7, 1:7))
for fig_id in 1:4
    f = Figure(size=(2000, 2000), fonts = (; regular = "Roboto"))
    for (i, ci) in enumerate(cidx)
        sp_idx = ((fig_id-1) * prod(size(cidx))) + i + 1
        if sp_idx <= length(spnames)
            plot_baseline_sdm(f, (ci[2], ci[1]), sdms[spnames[sp_idx]], title=spnames[sp_idx])
        end
    end 
    f
    save(joinpath(outdir, "sdm$fig_id.png"), f)
end
f
