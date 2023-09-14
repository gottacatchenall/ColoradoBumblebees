using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie, GeoMakie, GeoJSON, ColorSchemes

projected_overlap_dirs = [joinpath(artifactdir(), "projected_overlap", x) for x in readdir(joinpath(artifactdir(), "projected_overlap"))]

proj_overlaps = ColoradoBumblebees.load.(projected_overlap_dirs)



# List of figures:

# Main text:
# ----------------------------------------------------------------------
# - F1 Baseline interaction-richess/sdm-uncertainty bivariate map
# - F2 Each bee species overlap distribution across plants taken as a slice from SSP370, 2050s
# - F3 The mean overlap relative to baseline across all species pairs over each decade, with each SSP in colors


function standardize_layer(layer)
    z = StatsBase.fit(ZScoreTransform, Float32.(layer.grid))
    g = StatsBase.transform(z, Float32.(layer.grid)) 
    newlayer = similar(layer)
    newlayer.grid .= g
    z, newlayer
end

function get_quantiles(z, layer)
    quantile(StatsBase.transform(z, Float32.(layer.grid)), [0.25, 0.5, 0.75])
end

function quantize_layer(layer, quantiles)
    newlayer = similar(layer)
    for i in eachindex(layer.grid)
        tmp = findfirst(x-> layer.grid[i] < x, quantiles)
        newlayer.grid[i] = !isnothing(tmp) ? tmp : length(quantiles)+1
    end
    newlayer
end

function make_quantized_layers(baseline, future)
    int_rich, uncert = baseline.interaction_richness, baseline.sdm_uncertainty
    (z_rich, std_rich), (z_uncert, std_uncert) = standardize_layer(int_rich), standardize_layer(uncert)

    quant_richness, quant_uncert = get_quantiles(z_rich, int_rich), get_quantiles(z_uncert, uncert)
    richness_quantized, uncert_quantized = quantize_layer(std_rich, quant_richness), quantize_layer(std_uncert, quant_uncert)
end

#
# F1 - Baseline interaction-richess/sdm-uncertainty bivariate map
#
function baseline_bivar(baseline; num_bins = 4)
    baseline_richness, baseline_uncert = baseline.interaction_richness, baseline.sdm_uncertainty
    
    richness_quantized, uncert_quantized = make_quantized_layers(baseline, nothing)
    
    encoding = ([num_bins*j+i for i in 1:num_bins, j in 0:(num_bins-1)])

    bivar_layer = similar(baseline_richness)
    for I in eachindex(baseline_uncert.grid)
        bivar_layer.grid[I] = encoding[Int32(uncert_quantized.grid[I]), Int32(richness_quantized.grid[I])]
    end
    bivar_layer
end

function load_geojsons()
    counties = GeoJSON.read(read(datadir("public", "geojsons", "usa", "counties.json"), String))
    state = GeoJSON.read(read(datadir("public", "geojsons", "usa", "states.json"), String))
    state, counties
end

function setup_axis(fig_slice, layer)
    long = EXTENT[:left], EXTENT[:right]
    lat = EXTENT[:bottom], EXTENT[:top]
    ga = GeoAxis(
            fig_slice;
            titlealign=:left,
            yticklabelsvisible=false,
            xticklabelsvisible=false,
            xgridvisible=false,
            ygridvisible=false,
            yticks=34:2:44,
            xticks=-110:2:-102,
            dest="+proj=ortho +lon_0=-105 +lat_0=30",
            lonlims=(EXTENT[:left], EXTENT[:right]),
            latlims= (EXTENT[:bottom], EXTENT[:top]),
    )
    l = Matrix{Float32}(layer.grid)

    state, counties = load_geojsons()

    lats = collect(LinRange(lat[1], lat[2], size(l, 1)))
    longs = collect(LinRange(long[1], long[2], size(l, 2)))
    blue_red_4 = vec(reshape([color(x) for x in ["#d3d3d3", "#cda2a2", "#c66c6e", "#be2427", "#adc3cb", "#a8959c", "#a2646a", "#9c2226", "#87b2c3", "#838995", "#7f5b66", "#791f24", "#5fa1ba", "#5c7c8e", "#595361", "#551c22"]], (4,4)))

    hm = heatmap!(ga, longs, lats, l', colormap=blue_red_4)
    poly!(
        ga,
        counties;
        strokecolor=:white,
        strokewidth=1.5,
        color=(:grey, 0),
        shading=false,
    )
    
    poly!(ga, state; strokecolor=:white, strokewidth=2, color=(:blue, 0), shading=false)

end 


begin 
    layer = baseline_bivar(proj_overlaps[1])
    f = Figure(resolution=(2000,3000))
    setup_axis(f[1,1], layer)
    f
end

save(plotsdir("bivar_attempt.svg"), f)

#
# - F2 Each bee species overlap distribution across plants taken as a slice from SSP245, 2050s
#

function bee_slice(proj_overlap)

    bee_species = unique(proj_overlap.cooccurrence_dataframe.bee)

    category_label = String[]
    data_vec = Float32[]

    for sp_name in bee_species
        this_df = filter(x->x.bee == sp_name, proj_overlap.cooccurrence_dataframe)
        for r in eachrow(this_df)
            push!(category_label, r.bee)
            push!(data_vec, r.mean_cooccurrence)
        end
    end

    rainclouds(category_label, data_vec)
end

bee_slice(proj_overlaps[13])


#
# - F3 The mean overlap relative to baseline across all species pairs over each decade, with each SSP in colors
#

futures = proj_overlaps[2:end]
 

_timespan(::ProjectedOverlap{T,S}) where {T,S} = T
_scenario(::ProjectedOverlap{T,S}) where {T,S} = S

xvals = Dict([t => i for (i,t) in enumerate(TIMESPANS[2:end])])

begin 
fig = Figure(resolution=(2000, 1000))
axes = [Axis(
        fig[1,x],
        limits=(0, 5, 0.5,1.5)
    ) for x in values(xvals)]
cols = Dict([
    SSP1_26=>colorant"#2e81e8",
    SSP2_45=>colorant"#555",
    SSP3_70=>colorant"#e63e3e",
])



for (i,f) in enumerate(reverse(futures))
    t, s = _timespan(f), _scenario(f)
    ax = axes[xvals[t]]
    y = f.cooccurrence_dataframe.mean_cooccurrence
    density!(ax, y, direction=:y, color=(cols[s], 0.4))
end

fig


fw = PolyElement(color = (cols[SSP1_26], 0.7),)
pp = PolyElement(color = (cols[SSP2_45], 0.7),)
hp = PolyElement(color = (cols[SSP3_70], 0.7),)

Legend(fig[:,0], width=220, [fw,pp,hp], ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
fig 

end



# Supplement:
# ----------------------------------------------------------------------
# S?? - SDM Uncertainty: Each row: year range, Each column: SSP
# S?? - Interaction Richness: Each row: year range, Each column: SSP
# S?? - Bivariate: Each row : year, Each column, SSP 


#
# S?? - SDM Uncertainty: Each row: year range, Each column: SSP
#


function find_projection(proj_overlaps, scen, time)
    proj_overlaps[findall(isequal(time), _timespan.(proj_overlaps)) ∩ findall(isequal(scen), _scenario.(proj_overlaps))][1]
end 


begin
    fig = Figure(resolution=(1400, 2400))

    g = GridLayout(fig[1,1])

    for (i,t) in enumerate(TIMESPANS[2:end])
        for (j, ssp) in enumerate([SSP1_26, SSP2_45, SSP3_70])
            ax = Axis(g[i,j])
            
            proj = find_projection(proj_overlaps, ssp, t)
            heatmap!(ax, proj.interaction_richness)
        end
    end

    fig
end