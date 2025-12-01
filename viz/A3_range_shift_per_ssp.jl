"""
    This script makes maps of the area gained/lost over the century.
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

colormap = [
    colorant"#0573f9",
    colorant"#37a6f0",
    colorant"#64B9F1",
    colorant"#aed9f5",
    colorant"#555",
    colorant"#de3e3eff",
    colorant"#df5050ff",
    colorant"#e06060ff",
    colorant"#d87d7dff",
]

bgcolor = :grey98


ssps = ["SSP126", "SSP245", "SSP370"]

species = sort(collect(keys(sdms)))
outdir = joinpath("plots", "range_shifts")
mkpath(outdir)

num_cols, num_rows = 3, 7
num_figs = Int(ceil(length(species)/(num_cols*num_rows)))
species_cursor = 1

for figidx in 1:num_figs
    f = Figure(size=(2000, 2400), figure_padding=50)
    g = GridLayout(f[1,1])

    for i in CartesianIndices((1:num_cols, 1:num_rows))
        h = g[i[2], i[1], Makie.GridLayoutBase.Outer()] = GridLayout()
        box = Box(
            h[0:1,1:3], 
            color = :grey98,
            alignmode = Outside(-20, -20, -20, -20),
            strokecolor = :transparent,
            cornerradius = 15,
        )
        Makie.translate!(box.blockscene, 0, 0, -100)
        for i in 1:3
            ax = Axis(
                h[1,i], 
                title=ssps[i],
                titlealign=:left,
                titlefont=:regular,
                aspect=1
            )
            hidedecorations!(ax)
            sdm_ts = get_sdm_timeseries(sdms, species[species_cursor], ssps[i])
            gnl = compute_gains_and_losses(sdm_ts)

            #cm = [colormap[c] for c in unique(nodata(gnl, iszero))]
                
            
            #heatmap!(ax, nodata(gnl, !iszero), colormap=[bgcolor])
            heatmap!(ax, nodata(gnl, 0), colormap=colormap, colorrange=(1,9))
            lines!(ax, county_poly, color=:grey20, linewidth=0.15)
            lines!(ax, state_poly, color=:grey20, linewidth=0.7)
        end
        Label(h[0,:], text = species[species_cursor], font=:bold_italic, halign=:left, fontsize = 20)

        species_cursor += 1

        if species_cursor > length(sdms)
            break
        end

    end


    try # try to skip nonexistant rows for the final fig    
        for c in 1:num_cols-1
            colgap!(g, c, Relative(0.03))
        end
        for r in 1:num_rows-1
            rowgap!(g, r, Relative(0.03))
        end 
        catch 
    end 
     
    Legend(
        g[num_rows+1, :],
        [PolyElement(color=c, markersize=20, strokecolor=:grey50, strokewidth=1) for c in colormap],
        [   
        
            "Gained 2021-2040",
            "Gained 2041-2060",
            "Gained 2061-2080",
            "Gained 2081-2100",
            "Persists",
            "Lost 2021-2040", 
            "Lost 2041-2060",
            "Lost 2061-2080",
            "Lost 2081-2100",
        ],
        labelsize=28,
        orientation=:horizontal,
        framevisible=false,
        backgroundcolor = :grey94,
        nbanks=3,
    )
    #f
    save(joinpath(outdir, "$figidx.png"), f)
end 

