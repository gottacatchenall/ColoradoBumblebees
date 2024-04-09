using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie, GeoMakie
using GeoJSON
using ColorSchemes

data = load_data()

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=52)
set_theme!(fontsize_theme)


function load_geojsons()
    counties = GeoJSON.read(read(datadir("public", "geojsons", "usa", "counties.json"), String))
    state = GeoJSON.read(read(datadir("public", "geojsons", "usa", "states.json"), String))
    state, counties
end

function setup_axis(fig_slice, sdm, title)
    ga = GeoAxis(
            fig_slice;
            title = title, 
            titlealign=:left,
            yticklabelsvisible=false,
            xticklabelsvisible=false,
            yticks=34:2:44,
            xticks=-110:2:-102,
            dest="+proj=ortho +lon_0=-105 +lat_0=30",
            lonlims=(EXTENT[:left], EXTENT[:right]),
            latlims= (EXTENT[:bottom], EXTENT[:top]),
    )
end 

function plot_sdm(fig_slice, sdm::SpeciesDistribution)
    ax = setup_axis(fig_slice, sdm, sdm.species.name )
    l = Matrix{Float32}(sdm.probability.grid)

    thres = sdm.fit_stats["threshold"]
   # l[findall(x-> x < thres, l)] .= 0
    lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
    longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
    heatmap!(ax, longs, lats, ones(size(l')), colormap=[:grey80])

    hm = heatmap!(ax, longs, lats, l'; colormap=:magma, colorrange=(thres,1))
    
    state, counties = load_geojsons()

    poly!(ax, state, strokecolor=:white, strokewidth=5, color=(:blue, 0))
    poly!(ax, counties, strokecolor=:white, strokewidth=2, color=(:blue, 0))
end


species_name = "Bombus pensylvanicus"

species_name = "Aconitum columbianum"

data = load_data()
sp = occursin("Bombus", species_name) ? bee(data, species_name) : plant(data, species_name)



_baseline_sdm_cont = load_sdm(sp, baseline(), Baseline)
_baseline_sdm = load_sdm(sp, baseline(), Baseline)
_2050s = [load_sdm(sp, TIMESPANS[5], ssp) for ssp in [SSP1_26, SSP2_45, SSP3_70]]
_2090s = [load_sdm(sp, TIMESPANS[end], ssp) for ssp in [SSP1_26, SSP2_45, SSP3_70]]

τ = _baseline_sdm.fit_stats["threshold"]
for sdm in vcat(_baseline_sdm, _2050s..., _2090s...)
    neg_idx = findall(x -> x < τ, sdm.probability.grid)
    pos_idx = findall(x -> x >= τ, sdm.probability.grid)

    sdm.probability.grid[neg_idx] .= 0
    sdm.probability.grid[pos_idx] .= 1
end

function make_change_map(base, middle, last)
    _base = convert(Bool, base.probability)
    _middle = convert(Bool, middle.probability)
    _last = convert(Bool, last.probability)

    lost50s = _base.grid .&& .!_middle.grid
    lost90s = _base.grid .&& .!_last.grid
    shared = _base.grid .&& _middle.grid .&& _last.grid 
    gained50s = .!_base.grid .& _middle.grid
    gained90s = .!_base.grid .& _last.grid

    @info length.(findall.([lost50s, lost90s, shared, gained50s, gained90s]))

    change_map = similar(base.probability)
    change_map.grid .= NaN
    change_map.grid[findall(lost50s)] .= 1
    change_map.grid[findall(lost90s)] .= 2
    change_map.grid[findall(shared)] .= 3
    change_map.grid[findall(gained90s)] .= 4
    change_map.grid[findall(gained50s)] .= 5

    change_map.grid
end

fig =Figure(resolution=(3300, 2000))

main_fig = fig[1,1]
plot_sdm(main_fig, _baseline_sdm_cont)


g = GridLayout(fig[1,2])

cm = ColorScheme([
        colorant"#e63e3e",
        colorant"#f58282",
        colorant"#555",
        colorant"#9dc4f5",
        colorant"#2e81e8",
])


titles = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"]

for (i, fig_idx) in enumerate([(1,1), (1,2), (2,1)])
    ax = setup_axis(g[fig_idx...], _baseline_sdm, titles[i])
    l = make_change_map(_baseline_sdm, _2050s[i], _2090s[i])
    lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
    longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
    heatmap!(ax, longs, lats, ones(size(l')), colormap=[:grey80])
    hm = heatmap!(ax, longs, lats, l'; colormap=cm)
    state, counties = load_geojsons()


    poly!(ax, state, strokecolor=:white, strokewidth=3, color=(:blue, 0))
    poly!(ax, counties, strokecolor=:white, strokewidth=1, color=(:blue, 0))
end 


lost50s = PolyElement(color = (colorant"#e63e3e", 0.7), markersize = 55)
lost90s = PolyElement(color = (colorant"#f58282", 0.7), markersize = 55)
shared = PolyElement(color = (colorant"#555", 0.7), markersize = 55) 
gained50s = PolyElement(color = (colorant"#9dc4f5", 0.7), markersize = 55)
gained90s = PolyElement(color = (colorant"#2e81e8", 0.7), markersize = 55)

#=Legend(
    g[2,3], 
    [lost50s, lost90s, shared, gained50s, gained90s], 
    ["Lost by 2050s", "Lost by 2090s", "Persistent Range", "Gained by 2050s", "Gained by 2090s"], 
    framevisible = false)=#


current_figure()

save(plotsdir("F??_example_SDM projections.png"), fig)