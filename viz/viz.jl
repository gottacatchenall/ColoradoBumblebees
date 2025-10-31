using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG

# List of Figures:
# Species and Interaction Richness with elevation underneath
# Bivariate Interaction richness and uncertainty map
# β_wn map, bivariate w/ interaction richness 
# Most and least overlap among interacting pairs 
# ROCAUCs & TSS for each species 


# Fig Ideas:

# Correlation between overlap change and range area change
# Elevation displacement by species 
# 

include("io.jl")
include("networks.jl")
CairoMakie.activate!(; px_per_unit=3)

us_states = getpolygon(PolygonData(GADM, Countries), country="USA", level=1)
us_counties = getpolygon(PolygonData(GADM, Countries), country="USA", level=2)
states_to_include = ["Colorado", "Utah", "Wyoming", "Nebraska", "NewMexico", "Arizona", "Oklahoma", "Texas"]

bbox = (left=-109.7, right=-101.8, bottom=34.5,top=42.5)
bbox_poly = SDT.SimpleSDMPolygons._get_polygon_from_bbox(bbox)

elevation = SDMLayer(RasterData(WorldClim2, Elevation); resolution=0.5, bbox...)  



county_polys = intersect(vcat(
    [us_counties["Level One" => s] for s in states_to_include]...
), bbox_poly)
state_polys = intersect(FeatureCollection(
    [us_states[s] for s in states_to_include],
), bbox_poly)


sdms = load_sdms()


int_richness = get_interaction_richness(sdms; threshold=true) 

f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hidespines!(ax)
heatmap!(ax, int_richness)
lines!(ax, county_polys, color=:white, linewidth=0.15)
lines!(ax, state_polys, color=:white, linewidth=0.75)

ax2 = Axis(f[1,2], aspect=DataAspect())
hidespines!(ax2)
heatmap!(ax2, get_total_species_richness(sdms; threshold=true))
lines!(ax2, county_polys, color=:white, linewidth=0.15)
lines!(ax2, state_polys, color=:white, linewidth=0.75)

f




# --------- Begin Fit diagnostics ---------

species_dirs = [joinpath("artifacts", x) for x in readdir("artifacts")]

bee_color =("#5ca0f2ff", 0.4)
plant_color = ("#9d5bc5ff", 0.2)

X = [occursin("Bombus", x) ? 1 : 1.2 for x in species_dirs]
cols = [x == 1 ? bee_color : plant_color for x in X]

μ_roc = [open(joinpath("artifacts", x, "metrics.json"), "r") do f return JSON.parse(f)["rocauc"]["mean"] end
 for x in readdir("artifacts/")] #|> mean
σ_roc = [open(joinpath("artifacts", x, "metrics.json"), "r") do f return JSON.parse(f)["rocauc"]["std"] end
 for x in readdir("artifacts/")] 
μ_tss = [open(joinpath("artifacts", x, "metrics.json"), "r") do f return JSON.parse(f)["tss"]["mean"] end
 for x in readdir("artifacts/")] #|> mean
σ_tss = [open(joinpath("artifacts", x, "metrics.json"), "r") do f return JSON.parse(f)["tss"]["std"] end
 for x in readdir("artifacts/")] 




# Plot
begin 
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

f = Figure(size = (1200, 500))
ax = Axis(
    f[1,1];
    ylabel="ROC-AUC",
    shared_axis_args...
)
limits!(ax, 0.85, 1.3, 0.4, 1.)
hlines!(ax,[0.5], linestyle=:dash, color=:grey40, linewidth=3)
rc1 = rainclouds!(ax, X, μ_roc; shared_raincloud_args... )
ax2 = Axis(
    f[1,2];
    ylabel="TSS",
    shared_axis_args...
)
limits!(ax2, 0.85, 1.3, 0.4, 1.)
rc2 = rainclouds!(ax2, X, μ_tss; shared_raincloud_args...)
f
end 
save("plots/fit.png", f)

# --------- End Fit Diagnostics ---------

# --- Species Richness & Int Richness ----
species_richness = get_total_species_richness(sdms, threshold=true)
bee_richness = get_bee_richness(sdms, threshold=true)
plant_richness = get_plant_richness(sdms, threshold=true)

begin 
f = Figure(size=(1200, 700))
g = GridLayout(f[1,1])
ax1 = Axis(
    g[1,1], 
    aspect=1,
    title="Plant Richness",
    titlealign=:left,
    titlesize=22,
)
ax2 = Axis(
    g[1,2], 
    aspect=1,
    title = "Bee Richness",
    titlealign=:left,
    titlesize=22,
)
hm_plant = heatmap!(ax1, plant_richness, colormap=:magma)
lines!(ax1, county_polys, color=:white, linewidth=0.3)
lines!(ax1, state_polys, color=:white, linewidth=1.25)
Colorbar(g[2,1], hm_plant, vertical=false, flipaxis = false, label="Richness", labelsize=20, ticks=0:20:140)

hm_bee = heatmap!(ax2, bee_richness, colormap=:magma)
lines!(ax2, county_polys, color=:white, linewidth=0.3)
lines!(ax2, state_polys, color=:white, linewidth=1.25)
Colorbar(g[2,2], hm_bee, vertical=false, flipaxis = false, label="Richness", labelsize=20, ticks=0:2:20)

rowsize!(g, 2, Relative(0.05))

f
end 
save("plots/species_richness.png", f)

# --- End Species Richness & Int Richness ----




# --- Compute Network Uniqueness ----
mw = get_metaweb()
local_pools = get_local_species_pools(sdms)
local_nets = get_local_networks(mw, local_pools)

βwn = zeros(size(local_nets))
for i in CartesianIndices(local_nets)
    if local_nets[i] != []
        βwn[i] = KGL08(betadiversity(βWN, local_nets[i], mw))
    end 
end

βwn_map = deepcopy(int_richness)
βwn_map.grid .= βwn
# --- End Compute Network Uniqueness ----



# --- Bivariate Helpers ----
function _palette(; low=colorant"#d3d3d3", high=colorant"#120fe3", breaks=3)
    breakpoints = LinRange(0.0, 1.0, breaks)
    return Makie.ColorSchemes.weighted_color_mean.(breakpoints, high, low)
end

function discretize(layer, n::Integer)
    return (x -> round(Int64, x)).(rescale(layer, 1, n))
end


function make_bivariate(P,U,; nbreaks=5, xlabel="", ylabel="")
    colormap1 = _palette(; high=colorant"#4dcfff", breaks=nbreaks)
    colormap2 = _palette(; high=colorant"#be2427", breaks=nbreaks)
    colormatrix = [ColorBlendModes.blend.(c1, c2; mode=BlendMultiply) for c1 in colormap1, c2 in colormap2]

    m1 = discretize(quantize(P), nbreaks)
    m2 = discretize(quantize(U), nbreaks)
    category = similar(m1)
    for i in eachindex(category)
        category[i] = LinearIndices(colormatrix)[m1[i], m2[i]]
    end

    f = Figure(size=(1200, 900))
    g = GridLayout(f[1,1])
    ax = Axis(g[1,1]; aspect=DataAspect(), xticklabelsize = 26, yticklabelsize=26) 
    heatmap!(ax, category, colormap=vec(colormatrix))

    ax_inset = Axis(
        g[1, 2],
        width=Relative(0.8),
        height=Relative(0.8),
        aspect=1,
        halign=0.4,
        valign=0.5,
        xticklabelsvisible=false,
        xticksvisible=false,
        yticklabelsvisible=false,
        yticksvisible=false,
        xlabel = xlabel,
        ylabel = ylabel,
        xlabelsize = 26,
        ylabelsize = 26,
    )
    heatmap!(ax_inset, colormatrix)
    #hidedecorations!(ax_inset)

    colsize!(g, 2, Relative(0.3))

    lines!(ax, county_polys, color=:white, linewidth=0.3)
    lines!(ax, state_polys, color=:white, linewidth=1.25)
    f
end


# --- End Bivariate Helpers ----


# ------- Int Richness vs. β-div -------
βwn_map_zerod = nodata(βwn_map, iszero)
ir_zerod = deepcopy(int_richness)
ir_zerod.indices .= βwn_map_zerod.indices

f = make_bivariate(
    ir_zerod, 
    βwn_map; 
    nbreaks=5,
    xlabel = "Interaction Richness",
    ylabel = "Network Uniqueness"
)
heatmap!(f.content[1],nodata(βwn_map, !iszero), colormap=[:grey10])
lines!(f.content[1], county_polys, color=:white, linewidth=0.15)
lines!(f.content[1], state_polys, color=:white, linewidth=0.75)
f
save("plots/int_richness_vs_uniqueness.png", f)


# ------- Int Richness vs. SDM Uncertainty -------
f = make_bivariate(
    Float64.(species_richness),
    int_richness;
    nbreaks=5,
    xlabel="Interaction Richness",
    ylabel="SDM Uncertainty",
)

f
save("plots/int_richness_vs_uncertainty.png", f)


# ------- Int Richness vs. Species Richness -------
f = make_bivariate(
    int_richness, 
    get_total_uncertainty(sdms); 
    nbreaks=5,
    xlabel="Interaction Richness",
    ylabel="SDM Uncertainty",
)


# ------- Highest vs. Lowest overlap -------

function compute_bee_overlap(sdms; threshold=true)
    mw = get_metaweb()
    total_size = length(findall(first(values(sdms))[:prediction].indices))

    bee_names = filter(x->occursin("Bombus", x), collect(keys(sdms)))
    overlaps = [[] for _ in bee_names]
    for (i,species) in enumerate(bee_names)
        base_sdm = Bool.(get_prediction(sdms[species], threshold))
        interacting_species = get_interacting_species(mw, species)
        for int_species in interacting_species
            if int_species ∈ keys(sdms)
                int_sdm = Bool.(get_prediction(sdms[int_species], threshold))
                overlap = length(findall(base_sdm & int_sdm)) / (length(findall(base_sdm)) + length(findall(int_sdm)))
                push!(overlaps[i], overlap)
            end 
        end
    end
    bee_names, overlaps
end 

bee_names, overlaps = compute_bee_overlap(sdms)


mygrays = Makie.ColorScheme([
    colorant"#e63e3e",
    colorant"#f58282",
    colorant"#555",
    colorant"#9dc4f5",
    colorant"#2e81e8",
])


sort_idx = sortperm(median.(overlaps))
sp = collect(keys(sdms))


# total range size of bee as color 

range_sizes = []
for b in bee_names
    range = Bool.(get_prediction(sdms[b], true))
    range_prop = length(findall(range))/length(findall(range.indices))
    push!(range_sizes, range_prop)
end


x,y, c = [], [], []
ct = 0
for (i,j) in enumerate(sort_idx)
    for z in overlaps[j] 
        push!(x,i)
        push!(y,z)
        push!(c, range_sizes[i])
    end
end

col
col = [(get(Makie.ColorSchemes.magma, i/(0.05+maximum(range_sizes))), 0.4) for i in c]



begin

f = Figure(size=(2800,1400), figure_padding=50)
ax = Axis(
    f[1,1],
    xticks=(1.15:16.15, bee_names[sort_idx]),
    xlabel="",
    ylabel = "Proportion of Shared Range",
    xticklabelsize=42,
    yticklabelsize=36,
    ylabelsize=50,
    xgridvisible=false,
    xticklabelrotation=π / 2,
)
limits!(ax, 0.5, 17, 0, 0.5)
rc = rainclouds!(
        ax, 
        x, 
        y, 
        plot_boxplots = true,
        clouds = nothing,
        cloud_width=0,
        markersize=28,
        side_nudge=0.3,
        jitter_width = 0.15,
        gap=-0.5,
        color=col
)
Colorbar(f[1,2], width=30, colorrange=(0,0.05+maximum(range_sizes)), ticks=(0:0.05:0.05+maximum(range_sizes)), ticklabelsize=26, label="Range Size", labelsize=40, colormap=:magma)
f
end
save("plots/shared_range.png", f)



function plot_gmm(fig, slice, results; title="")
    x = results["data"]["x"]
    y = results["data"]["y"]
    best_k = results["best_model"]["n_components"]
    components = results["best_model"]["components"]
    waics = results["model_comparison"]["waics"]
    

    begindoy = (Date(2025, 5, 1) - Date(2025, 1, 1)).value
    enddoy = (Date(2025, 10, 1) - Date(2025, 1, 1)).value

    row, col = slice

    isfirstcol = col == 1
    islastrow = row == 7

    ylabel = isfirstcol ? "Number of Observations" : ""

    # Create visualization
    #fig = Figure(size=(900, 500))
  
    months = 5:9
    doy_range_by_month = [((Date(2025, m, 1) - Date(2025, 1, 1)).value, (lastdayofmonth(Date(2025, m,1)) - Date(2025, 1, 1)).value) for m in months]
    xticks = ([median(doy_range_by_month[i]) for i in eachindex(months)], [monthname(i) for i in months])

    # Plot 2: Best fit with 95% CI
    ax = Axis(
        fig[slice...],
        title=title,
        titlealign=:left,
        titlefont=:bold_italic,
        aspect=1,
        xlabel="",
        ylabel=ylabel,
        xticks = xticks,
        xticksvisible = false,
        xticklabelsvisible = islastrow,
        xticklabelrotation=π/2,
        xgridvisible=false,
        ygridvisible=false,
        #yticksvisible = isfirstcol,
        #yticklabelsvisible = isfirstcol,
    )
    ymax = maximum(y) + (0.05*maximum(y))
    xlims!(begindoy, enddoy)
    ylims!(0, ymax)
    # Generate prediction curves with CI using posterior samples
    x_pred = range(begindoy, enddoy, length=365*2)
    

    cols = [ :grey98, :grey85,]
    for (i,doy_range) in enumerate(doy_range_by_month)
        poly!(ax, Point2f[(doy_range[1], 0), (doy_range[1], ymax), (doy_range[2]+1, ymax), (doy_range[2]+1, 0)], color = (cols[1 + i % 2], 0.2))
    end 

    # Extract posterior samples for all components
    n_posterior_samples = length(components[1]["mu_samples"])
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
    
    # Compute mean and 95% CI
    y_pred_mean = vec(mean(y_pred_samples, dims=1))
    y_pred_lower = [quantile(y_pred_samples[:, j], 0.025) for j in 1:length(x_pred)]
    y_pred_upper = [quantile(y_pred_samples[:, j], 0.975) for j in 1:length(x_pred)]
    

    

    # Plot 95% CI band
    band!(ax, x_pred, y_pred_lower, y_pred_upper, 
          color=(:dodgerblue, 0.4), 
          label="95% CI"
    )
    
    # Plot data
    scatter!(ax, x, y, label="Data", markersize=7, color=(:black, 0.5))
    
    # Plot individual components using posterior mean by computing mean of component curves
    colors = [:red, :green, :orange, :purple,]
    for (idx, comp) in enumerate(components)
        # Compute component curves for each posterior sample
        y_component_samples = zeros(n_posterior_samples, length(x_pred))
        
        for s in 1:n_posterior_samples
            μ_s = comp["mu_samples"][s]
            σ_s = comp["sigma_samples"][s]
            A_s = comp["amplitude_samples"][s]
            
            for (j, x_val) in enumerate(x_pred)
                y_component_samples[s, j] = A_s * exp(-(x_val - μ_s)^2 / (2 * σ_s^2))
            end
        end
        
        # Take mean across posterior samples
        y_component_mean = vec(mean(y_component_samples, dims=1))
        
        lines!(
            ax, 
            x_pred, 
            y_component_mean, 
            label="Component $(comp["component_id"])",
            linestyle=:dash,
            linewidth=2,
            color=(colors[mod1(idx, length(colors))], 0.6)
        )
    end
    
    # Plot mean prediction
    lines!(
        ax, 
        x_pred, 
        y_pred_mean, 
        label="Posterior Mean",
        linewidth=3,
        color=:dodgerblue
    )
    
    #axislegend(ax, position=:rt)
    return fig
end

spnames = readdir("artifacts")
paths = [joinpath("artifacts", sp, "phenology.json") for sp in spnames]

cidx = CartesianIndices((1:7, 1:7))

for fig_id in 1:4
#fig_id = 1
f = Figure(size=(2000, 2000), fonts = (; regular = "Roboto"))
for (i, ci) in enumerate(cidx)
    sp_idx = ((fig_id-1) * prod(size(cidx))) + i + 1
    if sp_idx <= length(spnames)
        res = load_gmm(paths[sp_idx])
        plot_gmm(f, (ci[2], ci[1]), res; title=spnames[sp_idx])
    end
end 
f
    save("plots/phen$fig_id.png", f)
end


sp_idx = 3
sp_name = "Bombus huntii"
res = load_gmm(joinpath("artifacts", sp_name, "phenology.json"))
f = Figure(fonts = (; regular = "Roboto"))
plot_gmm(f, (1,1), res; title=sp_name)
f 






