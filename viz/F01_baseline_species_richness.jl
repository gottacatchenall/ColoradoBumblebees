"""
    This script plots baseline species richness maps for bees and plants separately.
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
include("../src/networks.jl")

CairoMakie.activate!(; px_per_unit=3)
sdms = read_sdms("./artifacts")

state_poly, county_poly = get_polygons()

species_richness = get_total_species_richness(sdms)
bee_richness = get_bee_richness(sdms)
plant_richness = get_plant_richness(sdms)


usa = getpolygon(PolygonData(OpenStreetMap, Places), place="USA")
us_states = getpolygon(PolygonData(GADM, Countries), country="USA", level=1)
us_states = us_states - us_states["Alaska"] - us_states["Hawaii"] 

elev = SDMLayer(RasterData(WorldClim2, Elevation); left=-130, right=-60, bottom = 20.) 
mask!(elev, usa)
elev = trim(elev)
BOUNDING_BOX = (left=-109.7, right=-101.8, bottom=34.5, top=42.5)
bbox_poly = SimpleSDMPolygons._get_polygon_from_bbox(BOUNDING_BOX)


species_richness = Float32.(get_total_species_richness(sdms))
uncertainty = get_total_uncertainty(sdms)

nbreaks = 4

R = discretize(quantize(species_richness, nbreaks), nbreaks)
U = discretize(quantize(uncertainty, nbreaks), nbreaks)


bivar, colormatrix = make_bivariate(
    Float32.(R),
    Float32.(U),
    #high1=colorant"#21baf7",
    high2=colorant"#486bb8",
    high1=colorant"#61e3a6",
    nbreaks=nbreaks
)

bivar |> unique

# area in each colorbar 

areas = []
unique_vals = unique(bivar)

for i in unique_vals
    push!(areas, length(findall(isequal(i), bivar))/prod(size(bivar)))
end

sortidx = sortperm(areas)[5:end]

areas[sortperm(areas)]




function makelab!(ax, start, end1, end2, label; path=Ann.Paths.Corner(), style=Ann.Styles.LineArrow(), labelspace=:data, kwargs...)
    #sx = (start[1] * cos(θ) - start[2] * sin(θ), start[2] * cos(θ) + start[1] * sin(θ))
    sx = start[1], start[2]
    e1x = end1[1], end1[2]
    e2x = end2[1], end2[2]

    annotation!(ax, sx..., e1x...;
        text=label,
        path=path,
        style=style,
        labelspace=labelspace,
        kwargs...
    )
    annotation!(ax, sx..., e2x...;
        text=label,
        path=path,
        style=style,
        labelspace=labelspace,
        kwargs...
    )
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


function get_histogram_per_elevation(
    target, 
    elevation;
    num_bins = 50
)
    counts = zeros(Int, num_bins)

    target_vec = [target[i] for i in eachindex(target)]
    elevation_vec = [elevation[i] for i in eachindex(target)]


    elev_min, elev_max = extrema(elevation_vec)
    bin_edges = range(elev_min, elev_max, length=num_bins+1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:num_bins]
    target_sum = zeros(num_bins)

    
    for (elev, targ) in zip(elevation_vec, target_vec)
        bin_idx = min(searchsortedfirst(bin_edges, elev), num_bins)
        if bin_idx > 0
            target_sum[bin_idx] += targ
            counts[bin_idx] += 1
        end
    end

    target_avg = [counts[i] > 0 ? target_sum[i] / counts[i] : 0.0 for i in 1:num_bins]
    return extrema(elevation_vec), bin_centers, target_avg
end

elev = get_elevation(sdms)


num_bins = 250
_, bin_centers, avg_unc = get_histogram_per_elevation(rescale(uncertainty), elev; num_bins=num_bins)
_, bin_centers, avg_richness = get_histogram_per_elevation(species_richness, elev; num_bins=num_bins)

begin
f = Figure(size=(1500, 1000))
g = GridLayout(f[1,1])
ax = Axis(
    g[1,1],
    aspect=DataAspect()
)
heatmap!(ax, bivar, colormap=vec(colormatrix))

lines!(ax, county_poly, color=:grey10, linewidth=0.5)
lines!(ax, state_poly, color=:grey10, linewidth=1.5)

g_right = GridLayout(g[1,2])

#=
ax_usa = GeoAxis(
    g_right[2,1],
    dest="+proj=ortho +lon_0=-105 +lat_0=30"
)
hidedecorations!(ax_usa)
heatmap!(ax_usa, elev, colormap=[:grey92])
lines!(ax_usa, us_states, color=:grey20, linewidth=0.5)
poly!(ax_usa, bbox_poly, strokecolor=:grey30, strokewidth=3, color=(:pink, 0.4))
=# 

ax_barplot = Axis(
    g_right[3,1],
    aspect=1,
    ylabel = "Proportion of Area",
    rightspinevisible=false,
    topspinevisible=false,
    xgridvisible=false,
    ygridvisible=false,
    xticksvisible=false,
    xticklabelsvisible=false
)
ylims!(ax_barplot, 0, 0.25)
barplot!(ax_barplot, areas[sortidx], color=vec(colormatrix)[unique_vals[sortidx]])


ax_legend = Axis(
    g_right[2,1],
    aspect=1,
    #xticksvisible=false,
    #yticksvisible=false,
    #xticklabelsvisible=false,
    #yticklabelsvisible=false,
    width=320,
    alignmode=Outside()
)

hidedecorations!(ax_legend)
hidespines!(ax_legend)

# Timslop below

xp = LinRange(-1, 1, size(colormatrix, 1) + 1)
yp = LinRange(-1, 1, size(colormatrix, 2) + 1)

θ = π / 4
for i in axes(colormatrix, 1)
    xc = (xp[i], xp[i+1]) .+ (0.015, -0.015)
    for j in axes(colormatrix, 2)
        yc = (yp[j], yp[j+1]) .+ (0.015, -0.015)
        corners = [(xc[1], yc[1]), (xc[2], yc[1]), (xc[2], yc[2]), (xc[1], yc[2])]
        r_corners = [
            (c[1] * cos(θ) - c[2] * sin(θ), c[2] * cos(θ) + c[1] * sin(θ))
            for c in corners
        ]
        poly!(ax_legend, r_corners, color=colormatrix[i, j], strokecolor=:black, strokewidth=0.5)
    end
end

makelab!(ax_legend, (-1, 1.2), (-0.2, 1.25), (-1.2, 0.3), "High\nUncertainty"; justification=:center)
makelab!(ax_legend, (1, 1.2), (0.2, 1.25), (1.2, 0.3), "High Species\nRichness"; justification=:center)
makelab!(ax_legend, (1, -1.2), (1.2, -0.3), (0.15, -1.25), "Low\nUncertainty"; justification=:center,)
makelab!(ax_legend, (-1, -1.2), (-1.22, -0.2), (-0.15, -1.25), "Low Species\nRichness"; justification=:center)


# elev / unc bars
ax = Axis(g_right[1,1], 
    xaxisposition = :bottom, 
    aspect=1, 
    xlabel = "Mean Uncertainty",
    xticks=0:0.2:1,
    ylabel = "Elevation (m)",
    xgridvisible=false,
    ygridvisible=false,
    xlabelfont=:bold,
    xlabelcolor = colorant"#385faeff",
    xticklabelcolor = colorant"#385faeff",
)
limits!(ax, 0, 1, 1000, 4200)
barplot!(ax, bin_centers, avg_unc, color=("#3764c5ff", 0.8), direction=:x, gap=0)
ax2 = Axis(
    g_right[1,1],
    aspect=1, 
    xlabel = "Mean Species Richness", 
    xaxisposition = :top,
    xgridvisible=false,
    ygridvisible=false,
    xlabelcolor = colorant"#2e8c60ff",
    xlabelfont=:bold,
    xticklabelcolor = colorant"#2e8c60ff"
)
limits!(ax2, 0, 150, 1000, 4200)
barplot!(ax2, bin_centers, avg_richness, color=("#41a878", 0.55), direction=:x, gap=0)




colsize!(g, 1, Relative(0.7))
rowsize!(g_right, 2, Relative(0.5))
rowgap!(g_right, 1, Relative(0.0001))



f
end 

save(joinpath("plots", "richness_vs_uncertainty.png"), f)

