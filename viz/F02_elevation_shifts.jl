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

state_poly, county_poly = get_polygons()
sdms = read_sdms(joinpath("artifacts"))

ssps = ["SSP126", "SSP245", "SSP370"]
years = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

SSP_LABELS = Dict(
    "SSP126" => "SSP 1-2.6",
    "SSP245" => "SSP 2-2.5",
    "SSP370" => "SSP 3-7.0"
)


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

function get_future_elevation_range(sdms, species, ssp, year)
    elevation.grid[findall(!iszero, sdms[species][:future][ssp][year][:range])]
end



# ===============================================================
# Get median elevation shifts
# ===============================================================
elevation = get_elevation(sdms)



baseline_elevations = Dict([sp=>median(get_baseline_elevation_distribution(sdms, sp)) for sp in keys(sdms)])
shifts = Dict()
for ssp in ssps
    shifts[ssp] = Dict()
    for yr in years
        shifts[ssp][yr] = []

        future_elevations = Dict([sp=>get_future_elevation_range(sdms, sp, ssp, yr) for sp in keys(sdms)])
        for sp in keys(sdms)
            if !isempty(future_elevations[sp])
                push!(shifts[ssp][yr], median(future_elevations[sp]) - median(baseline_elevations[sp]))
            end 
        end 
    end 
end


ssp = "SSP245"
new_species = []
lost_species = []
for yr in years
    new = Float64.(similar(first(sdms)[2][:baseline][:range]))
    new.grid .= 0
    
    lost = Float64.(similar(first(sdms)[2][:baseline][:range]))
    lost.grid .= 0

    for (k,v) in sdms
        new += .!Bool.(v[:baseline][:range]) .& Bool.(v[:future][ssp][yr][:range])
        lost += Bool.(v[:baseline][:range]) .& .!Bool.(v[:future][ssp][yr][:range])
    end 
    push!(new_species, new)
    push!(lost_species, lost)
end


function quantize_with_reference(layers, reference; nbreaks=5)
    qs = quantile(values(nodata(reference, 0)), [i for i in LinRange(0,1,nbreaks+1)[2:end-1] ])
    quantized_layers = []
    for l in layers
        nd = nodata(l, 0)
        quantized = similar(nd)
        quantized.indices .= 1
        quantized.grid .= 0
        for i in eachindex(nd)
            firstval = findlast(x -> nd[i] >= x, qs)
            quantized.grid[i] = isnothing(firstval) ? 1 : firstval + 1
        end
        push!(quantized_layers, nodata(quantized, 0))
    end
    return quantized_layers, qs
end


function make_bivariate(QG, QL; nbreaks=5, high2=colorant"#759d77", high1=colorant"#6676a1")
    colormap1 = _palette(; high=high1, breaks=nbreaks)
    colormap2 = _palette(; high=high2, breaks=nbreaks)
    colormatrix = [ColorBlendModes.blend.(c1, c2; mode=BlendMultiply) for c1 in colormap1, c2 in colormap2]

    category = deepcopy(QG)
    for i in eachindex(QG)
        if !isnothing(QG[i]) && !isnothing(QL[i])
            category[i] = LinearIndices(colormatrix)[QG[i], QL[i]]
        else
            category.indices[i] = 0
        end
    end
    return category, colormatrix
end 

time_idx = 4
nbreaks = 4

ref = lost_species[end]
qs = quantile(values(nodata(ref, 0)), [i for i in LinRange(0,1,nbreaks+1)[2:end-1]])

quantized_loss, loss_qs = quantize_with_reference(lost_species, lost_species[end]; nbreaks=nbreaks)
quantized_gain, gain_qs = quantize_with_reference(new_species, new_species[end]; nbreaks=nbreaks)


quantized_loss[4] |> unique

bivar_axis_settings = (;
    xgridvisible=false,
    ygridvisible=false,
    titlealign=:left,
    titlesize=20,
)


TIME_COLORS = [
    (colorant"#3a7281ff", 0.3),
    (colorant"#5199c9", 0.3),
    (colorant"#3c58a3ff", 0.3),
    (colorant"#6a45acff", 0.3)
]
offset_factor = 200

begin 
    f = Figure(size=(1350,800))
    g = GridLayout(f[1,1])
    g1 = GridLayout(g[1,:])


    g2 = GridLayout(g[1,2])
  
    for (i, cidx) in enumerate([(1,1), (1,2), (2,1), (2,2)])
        ax = Axis(
            g2[cidx...];
            aspect=DataAspect(),
            title = years[i],
            bivar_axis_settings...
        )
        category, colormatrix = make_bivariate(
            Int.(quantized_gain[i]), Int.(quantized_loss[i]),
            high1=colorant"#3ec0fd",
            high2=colorant"#f65c5c",
            nbreaks=nbreaks
        )
        heatmap!(ax, category, colormap=vec(colormatrix))
        lines!(ax, county_poly, color=:grey30, linewidth=0.3)
        lines!(ax, state_poly, color=:grey30, linewidth=1.25)
    end 


    # === begin bivariate legend 
    ax_legend = Axis(
        g1[2,2],
        aspect=1,
        #xticksvisible=false,
        #yticksvisible=false,
        #xticklabelsvisible=false,
        #yticklabelsvisible=false,
        width=350,
        alignmode=Outside()
    )

    hidedecorations!(ax_legend)
    hidespines!(ax_legend)
    limits!(ax_legend, -1.7, 1.7, -1.7, 1.7)
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

    makelab!(ax_legend, (-1.2, 1.2), (-1.3, 0.2), (-0.3, 1.25), "Many\nSpecies\n Lost"; justification=:center)
    makelab!(ax_legend, (1.2, 1.2), (1.3, 0.2), (0.3, 1.25), "Many\nSpecies\nGained"; justification=:center)
    makelab!(ax_legend, (-1.2, -1.2), (-1.3, -0.2), (-0.3, -1.2), "Few\nSpecies\nGained"; justification=:center,)
    makelab!(ax_legend, (1.2, -1.2), (1.3, -0.2), (0.3, -1.2), "Few\nSpecies\nLost"; justification=:center)

    Legend(
        g1[2,1], 
        [MarkerElement(color = c, marker = :rect, markersize=40) for c in TIME_COLORS ], 
        years, 
        "Timespan", 
        position = :rt,
        width=180,
        rowgap = 10,
        backgroundcolor=:grey97,
        framewidth=0,
        padding=15,
    )

    # === end legend 



    ax = Axis(
        g1[1,:], 
        ylabel="Shift of Median Elevation (m)",
        yticks = -200:200:1200,
        xgridvisible=false,
        xticksvisible=false,
        xticklabelsvisible=false,
        ygridvisible=false,
        yticklabelfont=:bold,
        yticklabelsize=16,
        ylabelsize=18,
        topspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
    )
    limits!(ax, 0.0048, 0.02, -250, 1000, )
    for i in eachindex(ssps)
        for j in 1:4
            shift_list = shifts[ssps[i]][years[j]]
            density!(ax, shift_list, offset=i/offset_factor - (0.03/offset_factor), color=TIME_COLORS[j], direction=:y
            )
            vlines!(ax, [i/offset_factor - (0.03/offset_factor)], color=:grey60)
        end 
    end
    hlines!(ax, [0], linestyle=:dash, color=:grey70, linewidth=3)
    text!(ax, 0.0175, -70, text="No Change", color=:grey50, fontsize=16)


    text!(ax, 0.0055, 900, text="SSP 1-2.6", font=:bold, color=:grey30, fontsize=23)
    text!(ax, 0.0105, 900, text="SSP 2-4.5", font=:bold, color=:grey30, fontsize=23)
    text!(ax, 0.0155, 900, text="SSP 3-7.0", font=:bold, color=:grey30, fontsize=23)




    colsize!(g1, 1, Relative(0.3))
    colgap!(g1, 1, Relative(0.01))

    rowsize!(g1, 1, Relative(0.6))
    colsize!(g, 1, Relative(0.4))

    f
end 

save(joinpath("plots", "range_shift.png"), f)
