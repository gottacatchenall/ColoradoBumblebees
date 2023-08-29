using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie, GeoMakie
using GeoJSON

data = load_data()

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=32)
set_theme!(fontsize_theme)

function load_geojsons()
    counties = GeoJSON.read(read(datadir("public", "geojsons", "usa", "counties.json"), String))
    state = GeoJSON.read(read(datadir("public", "geojsons", "usa", "states.json"), String))
    state, counties
end

function setup_axis(fig_slice, sdm)
    ga = GeoAxis(
            fig_slice;
            subtitle="$(sdm.species.name)",
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
    ax = setup_axis(fig_slice, sdm)
    thres = sdm.fit_stats["threshold"]
    
    
    l = Matrix{Float32}(sdm.probability.grid)

    l[findall(x-> x < thres, l)] .= 0
    lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
    longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
    hm = heatmap!(ax, longs, lats, l'; colormap=:magma, colorrange=(0.4,1))
    
    state, counties = load_geojsons()

    poly!(ax, state, strokecolor=:white, strokewidth=5, color=(:blue, 0))
    poly!(ax, counties, strokecolor=:white, strokewidth=2, color=(:blue, 0))
end

function plot_uncertainty(sdm::SpeciesDistribution)
    fig = Figure(;)
    ax = setup_axis(fig, sdm)
    l = Matrix{Float32}(sdm.uncertainty.grid)
    lats = collect(LinRange(EXTENT[:bottom], EXTENT[:top], size(l, 1)))
    longs = collect(LinRange(EXTENT[:left], EXTENT[:right], size(l, 2)))
    hm = heatmap!(ax, longs, lats, l'; colormap=:magma, colorrange=(0,1))
    fig
end


# Each fig will be 4 rows, 3 cols?

allspecies = vcat(bees(data), plants(data))
allspecies = allspecies[sortperm(map(x->x.name, allspecies))]

n_batches = (length(allspecies) ÷ 16) + 1 

startidx = 1
for b in 1:n_batches
    startidx +=1 
    endidx = min(startidx + 16, length(allspecies) )

    these_species = allspecies[startidx:endidx]

    f = Figure(resolution=(2400, 3600))
    done = false

    ct = 1
    for i in 1:4  
        for j in 1:4
            @info these_species[ct]
            this_sdm = load_sdm(these_species[ct], baseline(), Baseline)
            plot_sdm(f[i,j], this_sdm)
            if these_species[ct] == these_species[end]
                done = true
                break
            end    
            ct += 1
            startidx += 1
        end

        if done 
            break
        end
    end 

    Colorbar(f[2:3,5], limits=(0.4,1), size=30, ticks=0.4:0.1:1, colormap = :magma, label="Probability")

    save(plotsdir("sdm_set$b.png"), f)
end 



current_figure()
save(plotsdir("sdm_test.png"), f)


#=
dodec = load_sdm(plant(load_data(), "Ericameria parryi"), baseline(), Baseline)
plot_sdm(Figure()[1,1], dodec)
current_figure()




these_species = allspecies[end-3:end]

f = Figure(resolution=(2400, 1200))
done = false

ct = 1
for i in 1:4  
        this_sdm = load_sdm(these_species[i], baseline(), Baseline)
        plot_sdm(f[1,i], this_sdm)
end 

Colorbar(f[:,5], limits=(0.4,1), size=30, ticks=0.4:0.1:1, colormap = :magma, label="Probability")
current_figure()
save(plotsdir("sdm_set11.png"), f)=#