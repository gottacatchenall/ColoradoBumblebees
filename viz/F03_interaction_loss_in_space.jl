"""
    Makes a map of the total number of lost interactions across space under each SSP,
    and a density of the median number of lost interactions across elevations.
"""

using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates
using Statistics

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG

include(joinpath("..", "src", "io.jl"))
include("shared.jl")
include("../src/networks.jl")

CairoMakie.activate!(; px_per_unit=3)

state_poly, county_poly = get_polygons()
sdms = read_sdms("./artifacts")
metaweb = get_metaweb()

SSP_LABELS = Dict(
    "SSP126" => "SSP 1-2.6",
    "SSP245" => "SSP 2-4.5",
    "SSP370" => "SSP 3-7.0"
)

num_bins = 100

N_QUANTILES = 5

SCENARIO_COLORS = [
    colorant"#EE977A",
    colorant"#F5745A", 
    colorant"#E44F3E",
]


QUANTILE_COLORS = vcat([get(Makie.ColorSchemes.Purples, i) for i in 0.2:0.15:0.95][begin:end-1]...)
years = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]
ssps = ["SSP126", "SSP245", "SSP370"]


function get_lost_interaction_map(sdms, bee, plant, ssp, year)
    baseline_bee = sdms[bee][:baseline][:range]
    baseline_plant = sdms[plant][:baseline][:range]

    future_plant = sdms[plant][:future][ssp][year][:range]
    future_bee = sdms[bee][:future][ssp][year][:range]

    baseline_int = Bool.(baseline_bee .& baseline_plant)
    future_int = Bool.(future_bee .& future_plant)

    loss = baseline_int .& .!future_int
    return loss
end

function get_lost_interaction_count(sdms, metaweb, ssp, year)
    ints = SpeciesInteractionNetworks.interactions(metaweb)

    loss_count = Float32.(similar(first(sdms)[2][:baseline][:range]))
    loss_count.grid .= 0
    for (bee, plant, _) in ints
        if bee ∈ keys(sdms) && plant ∈ keys(sdms)
            loss_count += get_lost_interaction_map(sdms, bee, plant, ssp, year)
        end
    end 
    return loss_count
end

function get_elevation(sdms)
    elevation = SDMLayer(RasterData(WorldClim2, Elevation); resolution=0.5)  
    template_sdm = first(sdms)[2][:baseline][:range]
    E, N =  eastings(template_sdm), northings(template_sdm)
    cropped_elevation = copy(template_sdm)
    for (i,e) in enumerate(E), (j,n) in enumerate(N)
        cropped_elevation.grid[j,i] = elevation[e,n]
    end
    return cropped_elevation
end 

loss_counts_by_ssp = Dict(
    [
        ssp=>[get_lost_interaction_count(sdms, metaweb, ssp, yr) for yr in years]
        for ssp in ssps
    ]
) 




# ========================================================
# Compute Lost interactions per quantile
# ========================================================

qs = quantile(
    values(nodata(loss_counts_by_ssp["SSP370"][end], 0)),
    [0.2i for i in 1:4]
)
quantized_lost = []

for ssp in ssps
    nd = nodata(loss_counts_by_ssp[ssp][end], 0)

    quantized = similar(nd)
    quantized.grid .= 0
    for i in eachindex(nd)
        firstval = findlast(x -> nd[i] >= x, qs)
        quantized.grid[i] = isnothing(firstval) ? 1 : firstval+1
    end
    push!(quantized_lost, quantized)
end 


function get_average_lost_interactions_per_elevation(
    loss, 
    elevation;
    num_bins = 50
)
    zerod_loss = nodata(loss, 0) 
    counts = zeros(Int, num_bins)

    interactions_vec = [zerod_loss[i] for i in eachindex(zerod_loss)]
    elevation_vec = [elevation[i] for i in eachindex(zerod_loss)]


    elev_min, elev_max = extrema(elevation_vec)
    bin_edges = range(elev_min, elev_max, length=num_bins+1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:num_bins]
    interactions_sum = zeros(num_bins)

    
    for (elev, interactions) in zip(elevation_vec, interactions_vec)
        bin_idx = min(searchsortedfirst(bin_edges, elev), num_bins)
        if bin_idx > 0
            interactions_sum[bin_idx] += interactions
            counts[bin_idx] += 1
        end
    end

    interactions_avg = [counts[i] > 0 ? interactions_sum[i] / counts[i] : 0.0 for i in 1:num_bins]
    return extrema(elevation_vec), bin_centers, interactions_avg
end


function plot_quantized_loss!(
    g,
    pos, 
    loss_map;
    title = "",
    halign=0.94
)
    ax = Axis(
        g[pos...],
        aspect=DataAspect(),
        title=title,
        xgridvisible=false,
        ygridvisible=false,
        backgroundcolor=:white,
        titlealign=:left,
        titlesize=18,
    )
    heatmap!(ax, loss_map, colormap=QUANTILE_COLORS, colorrange=(1, N_QUANTILES))
    lines!(ax, county_poly, color=:grey40, linewidth=0.3)
    lines!(ax, state_poly, color=:grey40, linewidth=1.25)

    inset_ax = Axis(
        g[pos...],
        width = Relative(0.2),
        height = Relative(0.2),
        halign = halign,
        valign = 0.,
        leftspinecolor=:grey50,
        rightspinecolor=:grey50,
        topspinecolor=:grey50,
        bottomspinecolor=:grey50,
        spinewidth=2,
        bottomspinevisible=false,
        rightspinevisible=false
    )
    areas = [length(findall(isequal(i), loss_map)) for i in 1:N_QUANTILES]
    limits!(inset_ax, 0.3, 5.8, 0, 1.05maximum(areas))
    poly!(inset_ax, [(0.3,0), (5.8,0), (5.8,1.05maximum(areas)), (0.3,1.05maximum(areas))], color=:white)
    barplot!(
        areas,
        color=QUANTILE_COLORS,
        strokewidth=1
    )
    hidedecorations!(inset_ax)
end


begin 
    f = Figure(size=(900, 900))
    g = GridLayout(f[1,1])

    plot_quantized_loss!(g, (1,1), quantized_lost[1], title="SSP 1-2.6")
    plot_quantized_loss!(g, (1,2), quantized_lost[2], title="SSP 2-4.5")

    bottom_grid = GridLayout(g[2,:])
    plot_quantized_loss!(bottom_grid, (1,1), quantized_lost[3], title="SSP 3-7.0", halign=0.96)


    Colorbar(
        bottom_grid[1,2],
        colormap=colormap = cgrad(QUANTILE_COLORS, N_QUANTILES, categorical = true),
        ticks = (0.5:4.5, string.(Int.(vcat(1,qs...)))),
        colorrange=(0.5,5.5),
        width=30,
        label = "Total Interactions Lost",
        flipaxis=false,
        labelsize = 16
    )

    # ----- Average Loss across elevation -----
    ax = Axis(
        bottom_grid[1,3], 
        xlabel = "Average Lost Interactions", 
        ylabel="Elevation (m)",
        xgridvisible = false,
        ygridvisible= false,
    )

    for i in 3:-1:1
        elev_extrema, bin_centers, avg_lost = get_average_lost_interactions_per_elevation(
                loss_counts_by_ssp[ssps[i]][4], 
                elev;
                num_bins=num_bins
        )
        x = avg_lost
        cutoff = findfirst(i -> i > x[end], x)

        band!(
            ax, 
            x[cutoff:end], 
            bin_centers[cutoff:end], 
            [bin_centers[end] for _ in bin_centers[cutoff:end]], 
            color=(SCENARIO_COLORS[i])
        )
        band!(
            ax,
            x[begin:cutoff], 
            bin_centers[begin:cutoff], 
            [maximum(bin_centers) for _ in bin_centers[begin:cutoff]], 
            color=(SCENARIO_COLORS[i]),  
            label = SSP_LABELS[ssps[i]],
        )
        lines!(
            ax, 
            x[cutoff:end], 
            bin_centers[cutoff:end],
            color=SCENARIO_COLORS[i],
            linewidth=3,
        )
        lines!(ax, x[1:cutoff], bin_centers[1:cutoff], color=SCENARIO_COLORS[i], linewidth=3)
    end 
    limits!(ax, 0, 420, 1000, 4200)

    colsize!(bottom_grid, 1, Relative(0.58))
    colgap!(bottom_grid, 2, Relative(0.05))

    axislegend(
        ax,
        position=:rb
    )

    f 
end 

save(joinpath("plots", "interaction_loss.png"), f)
