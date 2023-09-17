# List of figures:

# Main text:
# ----------------------------------------------------------------------
# - F004_baseline_bivariate.png Baseline interaction-richess/sdm-uncertainty bivariate map
# - F005_bee_ove rlap_midcentury.png Each bee species overlap distribution across plants taken as a slice from SSP370, 2050s
# - F006_overlap_over_time.png The mean overlap relative to baseline across all species pairs over each decade, with each SSP in colors

# Supplement:
# ----------------------------------------------------------------------
# S?? - SDM Uncertainty: Each row: year range, Each column: SSP
# S?? - Interaction Richness: Each row: year range, Each column: SSP
# S?? - Relative Interaction Richness
# S?? - Bivariate: Each row : year, Each column, SSP 
# S(??-??) - Bee overlap raincloud -- each row is a decade, 3 figs for each SSP



using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie, GeoMakie, GeoJSON, ColorSchemes
using HypothesisTests

projected_overlap_dirs = [joinpath(artifactdir(), "projected_overlap", x) for x in readdir(joinpath(artifactdir(), "projected_overlap"))]
proj_overlaps = ColoradoBumblebees.load.(projected_overlap_dirs)

best_rep = "RelativeAbundance_Structural_Environment"
best_fit_dir = joinpath(artifactdir(), "classification_fits", "multiple_representations", "brt", best_rep)

data = load_data()
bee_species, plants_species = bees(data), plants(data)
bee_species, plants_species = bee_species[sortperm([b.name for b in bee_species])], plants_species[sortperm([p.name for p in plants_species])]
binary_prediction, probability_prediction, empirical = ColoradoBumblebees.get_metaweb(best_fit_dir)
M = BipartiteNetwork( Matrix{Bool}(any.(binary_prediction .∪ empirical)), [string(b.name) for b in bee_species], [string(p.name) for p in plants_species],)

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=32)
set_theme!(fontsize_theme)


function standardize_layer(layer)
    z = StatsBase.fit(ZScoreTransform, Float32.(layer.grid))
    g = StatsBase.transform(z, Float32.(layer.grid)) 
    newlayer = similar(layer)
    newlayer.grid .= g
    z, newlayer
end

function get_quantiles(z, layer; num_bins=5)
    quantile(StatsBase.transform(z, Float32.(layer.grid)), [x for x in LinRange(0,1, num_bins+1)][2:end-1])         
    #quantile(StatsBase.transform(z, Float32.(layer.grid)), [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875])
    #quantile(StatsBase.transform(z, Float32.(layer.grid)), [0.2, 0.4,0.6,0.8])
end

function quantize_layer(layer, quantiles)
    newlayer = similar(layer)
    for i in eachindex(layer.grid)
        tmp = findfirst(x-> layer.grid[i] < x, quantiles)
        newlayer.grid[i] = !isnothing(tmp) ? tmp : length(quantiles)+1
    end
    newlayer 
end

function get_baseline_quantiles(baseline)
    int_rich, uncert = baseline.interaction_richness, baseline.sdm_uncertainty
    (z_rich, std_rich), (z_uncert, std_uncert) = standardize_layer(int_rich), standardize_layer(uncert)
end 

function make_quantized_baseline(baseline; num_bins = 5)
    (z_rich, std_rich), (z_uncert, std_uncert) = get_baseline_quantiles(baseline)
  
    quant_richness, quant_uncert = get_quantiles(z_rich, baseline.interaction_richness; num_bins=num_bins), get_quantiles(z_uncert, baseline.sdm_uncertainty; num_bins=num_bins)
    richness_quantized, uncert_quantized = quantize_layer(std_rich, quant_richness), quantize_layer(std_uncert, quant_uncert)
end

function make_quantized_future(baseline, future; num_bins=5)
    (z_rich, a), (z_uncert, _) = get_baseline_quantiles(baseline)
    quant_richness, quant_uncert = get_quantiles(z_rich, baseline.interaction_richness; num_bins=num_bins), get_quantiles(z_uncert, baseline.sdm_uncertainty; num_bins=num_bins)
    
    g = StatsBase.transform(z_rich, Float32.(future.interaction_richness.grid)) 
    std_rich = similar(a)
    std_rich.grid .= g

    g = StatsBase.transform(z_uncert, Float32.(future.sdm_uncertainty.grid)) 
    std_uncert = similar(a)
    std_uncert.grid .= g
    richness_quantized, uncert_quantized = quantize_layer(std_rich, quant_richness), quantize_layer(std_uncert, quant_uncert)
end


#
# F1 - Baseline interaction-richess/sdm-uncertainty bivariate map
#
function baseline_bivar(baseline; num_bins = 5)
    baseline_richness, baseline_uncert = baseline.interaction_richness, baseline.sdm_uncertainty
    
    richness_quantized, uncert_quantized = make_quantized_baseline(baseline;)
    
    @info unique(richness_quantized.grid), unique(uncert_quantized.grid)

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
    blue_red_8 = [color(x) for x in ["#d3d3d3", "#cfbcbd", "#c9a6a6", "#c59091", "#c07879", "#bb6061", "#b44446", "#ad2123", "#c2cacd", "#beb4b8", "#b99fa1", "#b48a8c", "#b07375", "#ab5b5e", "#a54144", "#9f1f22", "#b0c1c7", "#acacb3", "#a7989d", "#a38389", "#a06e72", "#9b575c", "#963e42", "#901e21", "#9eb8c1", "#9ba4ad", "#979198", "#937d84", "#90696f", "#8c5359", "#873b40", "#821d20", "#8cafbb", "#899ca8", "#858a94", "#827780", "#7f636b", "#7c4f56", "#78383e", "#731b20", "#7ba5b5", "#7893a2", "#75828f", "#72707c", "#705e68", "#6d4b54", "#69353c", "#651a1e", "#699cb0", "#668b9d", "#647b8b", "#616a78", "#5f5965", "#5c4651", "#59323a", "#56181e", "#5692a9", "#558398", "#527485", "#506474", "#4f5361", "#4c424e", "#4a2f38", "#47171c"]]
    blue_red_5 = vec(reshape([color(x) for x in ["#d3d3d3", "#cbacad", "#c28485", "#b9595b", "#ad2123", "#b4c3c9", "#ad9fa5", "#a57a7f", "#9e5257", "#941e22", "#95b3be", "#8f929c", "#897078", "#834c52", "#7a1c20", "#76a3b4", "#728594", "#6d6671", "#67454e", "#61191e", "#5692a9", "#53778b", "#4f5c6a", "#4c3e49", "#47171c"]], (5,5)))


    hm = heatmap!(ga, longs, lats, l', colormap=blue_red_5)
    poly!(
        ga,
        counties;
        strokecolor=:white,
        strokewidth=1.5,
        color=(:grey, 0),
        shading=false,
    )
    
    poly!(ga, state; strokecolor=:white, strokewidth=4, color=(:blue, 0), shading=false)

end 


begin 
    layer = baseline_bivar(proj_overlaps[1])
    f = Figure(resolution=(2000,3000))
    setup_axis(f[1,1], layer)
    f
end

save(plotsdir("F004_baseline_bivariate.svg"), f)
save(plotsdir("F004_baseline_bivariate.png"), f)

#
# - F2 Each bee species overlap distribution across plants taken as a slice from SSP245, 2050s
#

function bee_slice(baseline, proj_overlap)
    ov = proj_overlap.cooccurrence_dataframe.mean_cooccurrence ./ baseline.cooccurrence_dataframe.mean_cooccurrence
 
    I = [M[String(r.bee), String(r.plant)] for r in eachrow(proj_overlap.cooccurrence_dataframe)]

    relative_df = DataFrame(bee=proj_overlap.cooccurrence_dataframe.bee[I], plant=proj_overlap.cooccurrence_dataframe.plant[I], overlap=ov[I])



    bee_species = unique(proj_overlap.cooccurrence_dataframe.bee)

    sorting_ov = proj_overlaps[13].cooccurrence_dataframe.mean_cooccurrence ./ baseline.cooccurrence_dataframe.mean_cooccurrence
    sorting_df = DataFrame(bee=proj_overlap.cooccurrence_dataframe.bee[I], plant=proj_overlap.cooccurrence_dataframe.plant[I], overlap=sorting_ov[I])
    sorted_idx = sortperm([mean(filter(x->x.bee == sp_name, sorting_df).overlap) for sp_name in bee_species])
    
    
    category_label = String[]
    data_vec = Float32[]

    species_is_signif = []


    colors = []

    mygrays = ColorScheme([
        colorant"#e63e3e",
        colorant"#f58282",
        colorant"#555",
        colorant"#9dc4f5",
        colorant"#2e81e8",
    ])
   
    mns = [mean(filter(x->x.bee == sp_name, relative_df).overlap) for sp_name in bee_species[sorted_idx]]

    c = [(get(mygrays, 1-(1.3-x)/(1.3-0.75)), 0.4) for x in mns]
    


    pthresh = 1e-5

    for (i,sp_name) in enumerate(bee_species[sorted_idx])
        this_df = filter(x->x.bee == sp_name, relative_df)
        pval = pvalue(OneSampleTTest(this_df.overlap .- 1.0))
        pval < pthresh ? push!(species_is_signif, true) : push!(species_is_signif, false)
        
        for r in eachrow(this_df)
            push!(category_label, r.bee)
            push!(data_vec, r.overlap)

            push!(colors, c[i])
        end

    end

    category_label, data_vec, species_is_signif, colors
end

begin 
    fig = Figure(resolution=(2800,1400))
    g = GridLayout(fig[1,1])


    ax = Axis(
        g[1, 1];
        xticklabelrotation=π / 2,
        yticks=0:0.1:2.25,
        titlealign=:right,
        title="2050-2059",
        subtitle="SSP 2-4.5",
        titlesize=50,
        ylabel="Projected Spatial Overlap Relative to Baseline",
        xticklabelsize=30,
    )
    xlims!(ax, 0, 20)
    ylims!(ax, 0.65, 1.45)

    x,y, signif, cols= bee_slice(proj_overlaps[1], proj_overlaps[13])

    hlines!(ax, [1.0], linewidth=3, linestyle=:dash, color=:grey20)
    rainclouds!(
        ax, 
        x, 
        y, 
        plot_boxplots=false,
        cloud_width=0.5,
        markersize=28,
        side_nudge=0.2,
        jitter_width = 0.18,
        gap=-0.3,
        xticklabelcolor=:blue,
        color=cols,
        colormap=cols,
        colorrange=(0.4, 1.75),
    )


    ax2 = Axis(g[2,1], 
        topspinevisible=false, 
        leftspinevisible=false, 
        bottomspinevisible=false,
        rightspinevisible=false,
    )
    xlims!(ax2, 0,20)
    ylims!(ax2, 0,1)
    hidedecorations!(ax2)

    for (i,s) in enumerate(signif)
        s && text!(ax2, i-0.05, 0.1, text="*", fontsize=50)
    end

    rowsize!(g, 2, Relative(0.07))


    fig
end

save(plotsdir("F005_bee_overlap_midcentury.png"), fig)
save(plotsdir("F005_bee_overlap_midcentury.svg"), fig)

#
# - F3 The mean overlap relative to baseline across all species pairs over each decade, with each SSP in colors
#
baseline_po = proj_overlaps[1]
futures = proj_overlaps[2:end]

_timespan(::ProjectedOverlap{T,S}) where {T,S} = T
_scenario(::ProjectedOverlap{T,S}) where {T,S} = S

xvals = Dict([t => i for (i,t) in enumerate(TIMESPANS[2:end])])

begin 
    fig = Figure(resolution=(2000, 1000))
    axes = [Axis(
            fig[1,x],
            limits=(0, 5, 0.7,1.3)
        ) for x in values(xvals)]
    cols = Dict([
        SSP1_26=>colorant"#2e81e8",
        SSP2_45=>colorant"#555",
        SSP3_70=>colorant"#e63e3e",
    ])

    baseline = proj_overlaps[1]

    I = [M[String(r.bee), String(r.plant)] for r in eachrow(baseline.cooccurrence_dataframe)]


    for (i,f) in enumerate(reverse(futures))
        t, s = _timespan(f), _scenario(f)
        ax = axes[xvals[t]]
        y = f.cooccurrence_dataframe.mean_cooccurrence[I] ./ baseline.cooccurrence_dataframe.mean_cooccurrence[I]
        density!(ax, y, direction=:y, color=(cols[s], 0.4))
    end

    fig


    fw = PolyElement(color = (cols[SSP1_26], 0.7),)
    pp = PolyElement(color = (cols[SSP2_45], 0.7),)
    hp = PolyElement(color = (cols[SSP3_70], 0.7),)

    Legend(fig[:,0], width=220, [fw,pp,hp], ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
    fig 
end

save(plotsdir("F006_overlap_over_time.png"), fig)
save(plotsdir("F006_overlap_over_time.svg"), fig)


# SUPPLEMENT FIGS
grouped_fig_params = (;resolution=(2400,1400))
grouped_axis_params = ()

function find_projection(proj_overlaps, scen, time)
    proj_overlaps[findall(isequal(time), _timespan.(proj_overlaps)) ∩ findall(isequal(scen), _scenario.(proj_overlaps))][1]
end 

function load_geojsons()
    counties = GeoJSON.read(read(datadir("public", "geojsons", "usa", "counties.json"), String))
    state = GeoJSON.read(read(datadir("public", "geojsons", "usa", "states.json"), String))
    state, counties
end

# Total richness
begin
    colormax = maximum(vcat([[extrema(x.interaction_richness)...] for x in proj_overlaps]...))

    fig = Figure(; grouped_fig_params...)

    state, counties = load_geojsons()

    hms = []
    for (i,t) in enumerate(TIMESPANS[2:end])
        for (j, ssp) in enumerate([SSP1_26, SSP2_45, SSP3_70])
            title = j == 1 ? string(t)[1:4]*"s" : ""
            ga = GeoAxis(
                fig[j,i],
                title=title,
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
            proj = find_projection(proj_overlaps, ssp, t)
            l = proj.interaction_richness.grid
            lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
            longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
        
            push!(hms, heatmap!(ga, longs, lats, l', colorrange=(0,colormax)))

                

            poly!(ga, state, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
            poly!(ga, counties, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))

        end
    end
    
    for (i,str) in enumerate(["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
       ax = Axis(fig[i,0])
       hidedecorations!(ax)
       hidespines!(ax)
       limits!(ax, 0,1, 0,1)
       text!(ax, 0.1, 0.4, text=str, fontsize=40)

    end
    Colorbar(fig[2,end+1], hms[end], width=25 ,label="Interaction Richness")
    fig
end
save(plotsdir("S031_interaction_richness.png"), fig)


# Relative richness
begin

    fig = Figure(;grouped_fig_params...)

    g = GridLayout(fig[1,1])

    state, counties = load_geojsons()

    hms = []
    for (i,t) in enumerate(TIMESPANS[2:end])
        for (j, ssp) in enumerate([SSP1_26, SSP2_45, SSP3_70])
            title = j == 1 ? string(t)[1:4]*"s" : ""
            ga = GeoAxis(
                fig[j,i],
                title=title,
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
            proj = find_projection(proj_overlaps, ssp, t)
            l = proj.interaction_richness.grid
            lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
            longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
        
            a = proj.interaction_richness / baseline_po.interaction_richness
            l = a.grid
            push!(hms, heatmap!(ga, longs, lats, l', colorrange=(0,2.5)))

                

            poly!(ga, state, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
            poly!(ga, counties, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))

        end
    end
          
    for (i,str) in enumerate(["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
        ax = Axis(fig[i,0])
        hidedecorations!(ax)
        hidespines!(ax)
        limits!(ax, 0,1, 0,1)
        text!(ax, 0.1, 0.4, text=str, fontsize=40)
 
    end
    Colorbar(fig[2,end+1], hms[end], width=25 ,label="Interaction Richness Relative to Baseline")


    fig
end
save(plotsdir("S032_int_richness_relative_to_baseline.png"), fig)
 
# Uncert
begin
    colormax = maximum(vcat([[extrema(x.sdm_uncertainty)...] for x in proj_overlaps]...))

    fig = Figure(;grouped_fig_params...)
    state, counties = load_geojsons()

    g = GridLayout(fig[1,1])
    hms = []
    for (i,t) in enumerate(TIMESPANS[2:end])
        for (j, ssp) in enumerate([SSP1_26, SSP2_45, SSP3_70])
            title = j == 1 ? string(t)[1:4]*"s" : ""
            ga = GeoAxis(
                fig[j,i],
                title=title,
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
            proj = find_projection(proj_overlaps, ssp, t)
            l = proj.sdm_uncertainty.grid

            lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
            longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
        
            push!(hms, heatmap!(ga, longs, lats, l', colorrange=(0,colormax)))

            poly!(ga, state, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
            poly!(ga, counties, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
        end
    end
      
    for (i,str) in enumerate(["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
        ax = Axis(fig[i,0])
        hidedecorations!(ax)
        hidespines!(ax)
        limits!(ax, 0,1, 0,1)
        text!(ax, 0.1, 0.4, text=str, fontsize=40)
 
     end
    Colorbar(fig[2,end+1], hms[end], width=25 ,label="SDM Uncertainty")

    fig
end

save(plotsdir("S033_sdm_uncertainty.png"),fig)

#
# S?? - Bivariate: Each row : year, Each column, SSP 
#
function future_bivar(baseline, future; num_bins = 5)    
    richness_quantized, uncert_quantized = make_quantized_future(baseline, future; num_bins = 5)
    
    encoding = ([num_bins*j+i for i in 1:num_bins, j in 0:(num_bins-1)])
    bivar_layer = similar(richness_quantized)
    for I in eachindex(uncert_quantized.grid)
        bivar_layer.grid[I] = encoding[Int32(uncert_quantized.grid[I]), Int32(richness_quantized.grid[I])]
    end
    bivar_layer
end

begin
    blue_red_5 = vec(reshape([color(x) for x in ["#d3d3d3", "#cbacad", "#c28485", "#b9595b", "#ad2123", "#b4c3c9", "#ad9fa5", "#a57a7f", "#9e5257", "#941e22", "#95b3be", "#8f929c", "#897078", "#834c52", "#7a1c20", "#76a3b4", "#728594", "#6d6671", "#67454e", "#61191e", "#5692a9", "#53778b", "#4f5c6a", "#4c3e49", "#47171c"]], (5,5)))
    baseline_po = proj_overlaps[1]
    fig = Figure(;grouped_fig_params...)
    state, counties = load_geojsons()

    g = GridLayout(fig[1,1])
    hms = []
    for (i,t) in enumerate(TIMESPANS[2:end])
        for (j, ssp) in enumerate([SSP1_26, SSP2_45, SSP3_70])
            title = j == 1 ? string(t)[1:4]*"s" : ""
            ga = GeoAxis(
                fig[j,i],
                title=title,
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
            proj = find_projection(proj_overlaps, ssp, t)
            l = future_bivar(baseline_po, proj; num_bins=5).grid
            
            lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
            longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
        
            push!(hms, heatmap!(ga, longs, lats, l', colormap=blue_red_5))

            poly!(ga, state, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
            poly!(ga, counties, strokecolor=:white, strokewidth=0.75, color=(:blue, 0))
        end
    end
      
    for (i,str) in enumerate(["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"])
        ax = Axis(fig[i,0])
        hidedecorations!(ax)
        hidespines!(ax)
        limits!(ax, 0,1, 0,1)
        text!(ax, 0.1, 0.4, text=str, fontsize=40)
 
     end
    fig
end

save(plotsdir("S034_bivariate_overlap_futures.png"), fig)

# S?? - Bee overlap raincloud -- each row is a decade, 3 figs for each SSP

# 

begin

pltnames = Dict(
    SSP1_26=>"S035_bee_overlap_ssp126.png",
    SSP2_45=>"S036_bee_overlap_ssp245.png",
    SSP3_70=>"S037_bee_overlap_ssp370.png",
)

for ssp in [SSP1_26, SSP2_45, SSP3_70]
    fig = Figure(resolution=(1500, 2200))
    g = GridLayout(fig[1,1])


    
    for (i,t) in enumerate(TIMESPANS[2:end])
        ax=  Axis(
            g[i,1],
            xticklabelrotation=π / 2,
            yticks=0:0.1:2.25,
            titlealign=:right,
            title=split(string(t), "_")[1]*"s",
            titlesize=20,
            xticklabelsize=24,
            yticklabelsize=16,
            ylabelsize=20,
            ylabel ="Relative Spatial Overlap",
            xticklabelsvisible= t == TIMESPANS[end]
        )
        xlims!(ax, 0, 20)
        ylims!(ax, 0.65, 1.5)
        proj = find_projection(proj_overlaps, ssp, t)

        x,y, signif, cols = bee_slice(proj_overlaps[1], proj)


        
        hlines!(ax, [1.0], linewidth=2, linestyle=:dash, color=:grey20)
        rainclouds!(
            ax, 
            x, 
            y, 
            plot_boxplots=false,
            cloud_width=0.5,
            markersize=8,
            side_nudge=0.2,
            jitter_width = 0.14,
            gap=-0.3,
            xticklabelcolor=:blue,
            color=cols,
            colorrange=(0.4, 1.75),
        )
    end 

    fig
    save(plotsdir(pltnames[ssp]), fig)
end 
end 