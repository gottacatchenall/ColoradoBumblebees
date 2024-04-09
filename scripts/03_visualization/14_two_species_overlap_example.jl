using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie, GeoMakie
using GeoJSON
using ColorSchemes

data = load_data()

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=70)
set_theme!(fontsize_theme)

function load_future_binary_layers(species_name)
    data = load_data()
    sp = occursin("Bombus", species_name) ? bee(data, species_name) : plant(data, species_name)

    _baseline_sdm = load_sdm(sp, baseline(), Baseline)
    _2050s = [load_sdm(sp, TIMESPANS[5], ssp) for ssp in [SSP1_26, SSP2_45, SSP3_70]]
    _2090s = [load_sdm(sp, TIMESPANS[end], ssp) for ssp in [SSP1_26, SSP2_45, SSP3_70]]


    τ = _baseline_sdm.fit_stats["threshold"]
    sdms = vcat(_baseline_sdm, _2050s..., _2090s...)
    for sdm in sdms
        neg_idx = findall(x -> x < τ, sdm.probability.grid)
        pos_idx = findall(x -> x >= τ, sdm.probability.grid)

        sdm.probability.grid[neg_idx] .= 0
        sdm.probability.grid[pos_idx] .= 1
    end

    return [_baseline_sdm, _2050s, _2090s]
end 

function load_geojsons()
    counties = GeoJSON.read(read(datadir("public", "geojsons", "usa", "counties.json"), String))
    state = GeoJSON.read(read(datadir("public", "geojsons", "usa", "states.json"), String))
    state, counties
end

function setup_axis(fig_slice, sdm, title)
    ga = GeoAxis(
            fig_slice;
            subtitle=title,
            titlealign=:center,
            yticklabelsvisible=false,
            xticklabelsvisible=false,
            yticks=34:2:44,
            xticks=-110:2:-102,
            dest="+proj=ortho +lon_0=-105 +lat_0=30",
            lonlims=(EXTENT[:left], EXTENT[:right]),
            latlims= (EXTENT[:bottom], EXTENT[:top]),
    )
end 


function plot_overlap(fig_slice, bee_sdm::SpeciesDistribution, plant_sdm::SpeciesDistribution; title="")
    ax = setup_axis(fig_slice, bee_sdm, title)
    
    bee_binary = Matrix{Bool}(bee_sdm.probability.grid)
    plant_binary = Matrix{Bool}(plant_sdm.probability.grid)



    lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(bee_binary, 1)))
    longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(bee_binary, 2)))
    heatmap!(ax, longs, lats, ones(size(bee_binary')), colormap=[colorant"#eceff4"])

    overlap = bee_binary .& plant_binary
    bee_only = bee_binary .& .!plant_binary
    plant_only = .!bee_binary .& plant_binary
    neither = .!bee_binary .& .!plant_binary

    layer = zeros(size(bee_binary))
    
    layer[findall(plant_only)] .= 1
    layer[findall(bee_only)] .= 2
    layer[findall(overlap)] .= 3
    layer[findall(neither)] .= NaN

    #cm = [colorant"#3689d6", colorant"#55a176", colorant"#de3c54"]
    #cm = [colorant"#8fbcbb", colorant"#ebcb8b", colorant"#e0707b"]
    cm = [colorant"#ebcb8b", colorant"#5e81ac", colorant"#eb5967"]

    hm = heatmap!(ax, longs, lats, layer', colormap=cm)
    
    state, counties = load_geojsons()

    poly!(ax, state, strokecolor=:white, strokewidth=5, color=(:blue, 0))
    poly!(ax, counties, strokecolor=:white, strokewidth=2, color=(:blue, 0))
end


#bee_name = "Bombus pensylvanicus"
#plant_name = "Linaria vulgaris"

bee_name = "Bombus fervidus"
plant_name = "Achillea millefolium"

bee_layers = load_future_binary_layers(bee_name)
plant_layers = load_future_binary_layers(plant_name)


bee_baseline, plant_baseline = bee_layers[1], plant_layers[1]

fig = Figure(resolution=(3200, 2000))

g = GridLayout(fig[1,1])

plot_overlap(g[1,1], bee_baseline, plant_baseline; title="Baseline")
plot_overlap(g[1,2], bee_layers[2][2], plant_layers[2][2]; title="2050s")
plot_overlap(g[1,3], bee_layers[3][2], plant_layers[3][2]; title="2090s")

p = PolyElement(color = (colorant"#8fbcbb", 0.7), markersize = 55)
b = PolyElement(color = (colorant"#ebcb8b", 0.7), markersize = 55)
both = PolyElement(color = (colorant"#e0707b", 0.7), markersize = 55)

Legend(g[2,:], width=1200, [b,p,both], [plant_name, bee_name, "Overlap"], framevisible = false)

rowsize!(g, 2, Relative(0.3))
current_figure()


save(plotsdir("F??_example_spatial_overlap.png"), fig)