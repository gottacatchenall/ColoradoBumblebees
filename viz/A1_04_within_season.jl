"""
    This script plots the within season interaction richness quantiles across space for each month (left panels), and the overall phenology (right)
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
sdms = read_sdms(joinpath("artifacts"))

state_poly, county_poly = get_polygons()

function get_posterior_phenology(species::String; kwargs...)
    get_posterior_phenology(load_gmm(joinpath("./artifacts", species, "phenology.json")); kwargs...)
end

function get_posterior_phenology(results::Dict; normalized=true)
    components = results["best_model"]["components"]
    n_posterior_samples = length(components[1]["mu_samples"])

    begindoy = (Date(2025, 5, 1) - Date(2025, 1, 1)).value
    enddoy = (Date(2025, 10, 1) - Date(2025, 1, 1)).value
    x_pred = collect(begindoy:enddoy)
    y_pred_samples = zeros(n_posterior_samples, length(x_pred))

    # Generate predictions for each posterior sample
    for s in 1:n_posterior_samples
        for (j, x_val) in enumerate(x_pred)
            # Initialize to zero for this prediction point
            y_pred_samples[s, j] = 0.0
            
            # Sum over all components using sampled parameters
            for comp in components
                μ_s = comp["mu_samples"][s]
                σ_s = comp["sigma_samples"][s]
                A_s = comp["amplitude_samples"][s]
                y_pred_samples[s, j] += A_s * exp(-(x_val - μ_s)^2 / (2 * σ_s^2))
            end
        end
    end
    
    y_pred_mean = vec(mean(y_pred_samples, dims=1))
    y_pred_lower = [quantile(y_pred_samples[:, j], 0.025) for j in 1:length(x_pred)]
    y_pred_upper = [quantile(y_pred_samples[:, j], 0.975) for j in 1:length(x_pred)]
    
    M = maximum(y_pred_mean)

    if normalized
        y_pred_mean ./= M
        y_pred_lower ./= M
        y_pred_upper ./= M
    end

    return x_pred, y_pred_mean, y_pred_lower, y_pred_upper
end

function get_month_doys(month_num)
    ((Date(2025, month_num, 1) - Date(2025, 1, 1)).value, (lastdayofmonth(Date(2025, month_num, 1)) - Date(2025, 1, 1)).value)
end

function get_median_month_phenology_score(doys, phenology, month)
    startdoy, enddoy = get_month_doys(month)
    idx = findall(x-> x >= startdoy && x <= enddoy, doys)
    median(phenology[idx])
end

function get_phenology_score_by_month(species, month)
    res = load_gmm(joinpath("./artifacts", species, "phenology.json"))
    doys, phenology, _, _ = get_posterior_phenology(res)
    phen_score = get_median_month_phenology_score(doys, phenology, month)
end

function get_phenology_weighted_range(sdms, species, month)
    phen_score = get_phenology_score_by_month(species, month)
    phen_score .* (sdms[species][:prediction] .> sdms[species][:metrics]["threshold"]["mean"])
end 

function get_summed_phenology(species)
    doys, ytemplate, _, _ = get_posterior_phenology(load_gmm(joinpath("./artifacts", species[begin], "phenology.json")))
    y_median = zeros(length(ytemplate))
    y_lower = zeros(length(ytemplate))
    y_upper = zeros(length(ytemplate))

    for sp in species
        res = load_gmm(joinpath("./artifacts", sp, "phenology.json"))
        _, phenology, lower, upper = get_posterior_phenology(res)
        y_median .+= phenology
        y_lower .+= lower
        y_upper .+= upper
    end 
    return doys, y_median, y_lower, y_upper
end

function get_summed_plant_phenology(sdms)
    plant_species = [x for x in keys(sdms) if !occursin("Bombus", x)]
    get_summed_phenology(plant_species) 
end
function get_summed_bee_phenology(sdms)
    bee_species = [x for x in keys(sdms) if occursin("Bombus", x)]
    get_summed_phenology(bee_species) 
end


function get_interaction_richness_by_month(sdms, month)
    int_richness = Float64.(copy(first(values(sdms))[:baseline][:range]))
    int_richness.grid .= 0
    metaweb_df = get_metaweb_df()
    taxa_df = get_taxa_df("./data")
    filter!(x-> x.bee_name ∈ taxa_df.species_name && x.plant_name ∈ taxa_df.species_name, metaweb_df)

    
    ranges = Dict([s=>sdms[s][:baseline][:range] for s in keys(sdms)])
    phen_scores = Dict([s=>get_phenology_score_by_month(s, month) for s in keys(sdms)])

    for r in eachrow(metaweb_df)
        if r.interacts
            b_range = ranges[r.bee_name]
            p_range = ranges[r.plant_name]
            b_score = phen_scores[r.bee_name]
            p_score = phen_scores[r.plant_name]
            int_richness.grid[findall(isone, b_range .& p_range)] .+= b_score * p_score;
        end
    end    
    return int_richness
end



# --------------------------------------------------------------------------------
# Seasonal elevation shift
# --------------------------------------------------------------------------------

elevation = get_elevation(sdms)

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

function get_baseline_elevation_distribution(sdms, species)
    elevation.grid[findall(!iszero, sdms[species][:baseline][:range])]
end

function get_weighted_elevation_average(species, lower=0.25, upper=0.75)
    median_elevation_per_species = Dict([s=>median(get_baseline_elevation_distribution(sdms, s)) for s in species])
    _lower_elevation_per_species = Dict([s=>quantile(get_baseline_elevation_distribution(sdms, s), upper) for s in species])
    _upper_elevation_per_species = Dict([s=>quantile(get_baseline_elevation_distribution(sdms, s), lower) for s in species])

    x, y_template = get_posterior_phenology(first(species))

    phen_dict = Dict()
    for s in species
        x,y = get_posterior_phenology(s)
        phen_dict[s] = y
    end 

    weighted_elevation_avg = Float64[]
    weighted_elevation_lower = Float64[]
    weighted_elevation_upper = Float64[]

    for i in eachindex(y_template)
        comp = Dict([k=>v[i] for (k,v) in phen_dict])
        sum_comp = sum(values(comp))
        normalized_comp = Dict([k=>v/sum_comp for (k,v) in comp])

        push!(weighted_elevation_avg, sum([median_elevation_per_species[s] * normalized_comp[s] for s in species]))
        push!(weighted_elevation_upper, sum([_upper_elevation_per_species[s] * normalized_comp[s] for s in species]))
        push!(weighted_elevation_lower, sum([_lower_elevation_per_species[s] * normalized_comp[s] for s in species]))
    end

    return weighted_elevation_avg, weighted_elevation_lower, weighted_elevation_upper
end 


x, y_plant, lower_plant, upper_plant = get_summed_plant_phenology(sdms)
x, y_bee, lower_bee, upper_bee = get_summed_bee_phenology(sdms)


bee_elev_avg, bee_elev_low, bee_elev_up = get_weighted_elevation_average(filter(contains("Bombus"), keys(sdms)))
plant_elev_avg, plant_elev_low, plant_elev_up = get_weighted_elevation_average(filter(!contains("Bombus"), keys(sdms)))

int_richnesses = [get_interaction_richness_by_month(sdms, m) for m in 5:9]
maximum.(int_richnesses)

colormap = [
    colorant"#2e3440", 
    colorant"#4c566a",
    colorant"#54b2c7", 
    colorant"#e9764c", 
    colorant"#ebcb8b"
]

#plant_color = "#2c8f5eff"
#plant_label_color = "#2c8f5eff"

plant_color = "#5e9e9dff"
plant_label_color = "#4d8382ff"

#bee_color = "#ef9c5d"
#bee_label_color = "#e28035ff"

bee_color = "#5f95c7ff"
bee_label_color = "#5f95c7ff"




months = 5:9
doy_range_by_month = [((Date(2025, m, 1) - Date(2025, 1, 1)).value, (lastdayofmonth(Date(2025, m,1)) - Date(2025, 1, 1)).value) for m in months]
xticks = ([median(doy_range_by_month[i]) for i in eachindex(months)], [monthname(i) for i in months])

begindoy = 120
enddoy= 272

plant_ymax = 130
bee_ymax = 19.5

begin 
    f = Figure(size=(1500, 1300))
    g = GridLayout(f[1,1])
    
    g_int_richness = GridLayout(g[1,1])
    
    cidxs = [(1,1), (1,2), (2,1), (2,2), (3,2)]
    monthnames = [monthname(i) for i in months] 
    for (i, ir) in enumerate(int_richnesses)
        ax = Axis(
            g_int_richness[cidxs[i]...],
            aspect=1,
            title = monthnames[i],
            titlealign = :left,
            titlesize=24,
        )
        #heatmap!(ax, nodata(ir, x->!iszero(x)), colormap=[:grey98])
        heatmap!(ax, quantize(nodata(ir, 0), 5), colormap=:navia)
        #heatmap!(ax, quantize(ir, 5), colormap=:navia)
        lines!(ax, county_poly, color=:grey20, linewidth=0.15)
        lines!(ax, state_poly, color=:grey20, linewidth=0.7)
    end

    Legend(
        g_int_richness[3,1],
        [MarkerElement(color = c, marker = :rect, markersize = 30, strokecolor = :grey50, strokewidth = 1) for c in [get(Makie.ColorSchemes.navia, x) for x in 0:0.25:1]],
        ["0-20%", "20-40%", "40-60%", "60-80%", "80-100%"],
        "Interaction\nRichness\nQuantiles",
        framevisible = false,
        backgroundcolor = :grey94,
        labelsize = 25,
        titlesize= 26,
    )

    colsize!(g_int_richness, 1, Relative(0.5))
    #colsize!(g_int_richness, 2, Relative(0.33))


    # ----------------------------------------------------------------------
    # Phenology
    # ----------------------------------------------------------------------
    g_phenology = GridLayout(g[1,2])
    ax1 = Axis(
        g_phenology[1, 1], 
        aspect=1,
        xticks=xticks, 
        xticksvisible=false, 
        xgridvisible=false, 
        ygridvisible=true, 
        ygridstyle=:dash,
        ylabelfont=:bold,
        yticklabelcolor = plant_label_color, 
        ylabel="Plant Richness", 
        ylabelcolor=plant_label_color,
        yticks = 0:20:120,
        ylabelsize = 24,
        xticklabelsize = 18,
        yticklabelsize=18,
    )
    limits!(ax1, begindoy, enddoy, 0, plant_ymax)

    ax2 = Axis(
        title = "Species Richness",
        titlealign=:left,
        g_phenology[1, 1], 
        aspect=1,
        ygridvisible=false,
        yticklabelcolor = bee_label_color, 
        yaxisposition = :right, 
        ylabel="Bee Richness", 
        ylabelfont=:bold,
        ylabelcolor=bee_label_color,
        yticks = 0:3:18,
        titlesize = 25,
        ylabelsize = 24,
        yticklabelsize=18,
    )
    limits!(ax2, begindoy, enddoy, 0, bee_ymax)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    cols = [ :grey98, :grey80,]
    for (i,doy_range) in enumerate(doy_range_by_month)
        poly!(ax2, Point2f[(doy_range[1], 0), (doy_range[1], bee_ymax), (doy_range[2]+1, bee_ymax), (doy_range[2]+1, 0)], color = (cols[1 + i % 2], 0.2))
    end 
    lines!(ax1, x, y_plant, color = plant_color, linewidth=4)
    lines!(ax2, x, y_bee, color = bee_color, linewidth=4)
    band!(ax2, x, lower_bee, upper_bee, color=(bee_color, 0.4))
    band!(ax1, x, lower_plant, upper_plant, color=(plant_color, 0.4))
   
    colsize!(g, 1, Relative(0.65))
    colgap!(g, 1, Relative(0.03))
    

    elev_min,elev_max =2200, 3100
    ax_elevation = Axis(
        g_phenology[2, 1], 
        aspect=1,
        xticks=xticks, 
        xticksvisible=false, 
        xgridvisible=false, 
        ygridvisible=true, 
        ygridstyle=:dash,
        ylabelfont=:bold,
        title="Weighted Average Elevation",
        titlealign=:left,
        titlesize=25,
        ylabel="Elevation (m)", 
        ylabelsize = 24,
        xticklabelsize = 18,
        yticklabelsize=18,
    )
    limits!(ax_elevation, begindoy, enddoy, elev_min, elev_max)

    for (i,doy_range) in enumerate(doy_range_by_month)
        poly!(ax_elevation, Point2f[(doy_range[1], 0), (doy_range[1], elev_max), (doy_range[2]+1, elev_max), (doy_range[2]+1, 0)], color = (cols[1 + i % 2], 0.2))
    end 
    lines!(ax_elevation, x, plant_elev_avg, color=plant_color, linewidth=5)
    band!(ax_elevation, x, plant_elev_low, plant_elev_up, color=(plant_color, 0.4))
    lines!(ax_elevation, x, bee_elev_avg, color=bee_color, linewidth=5)
    band!(ax_elevation, x, bee_elev_low, bee_elev_up, color=(bee_color, 0.3))

    f
    #save(joinpath("plots", "within_season.png"), f)
end
