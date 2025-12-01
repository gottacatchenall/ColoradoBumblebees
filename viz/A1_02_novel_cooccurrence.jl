
# ----------------------------------------------------------------
# Novel Cooccurrence vs. Interactions Lost 
#
# ----------------------------------------------------------------


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

state_poly, county_poly = get_polygons()

sdms = read_sdms("./artifacts")

metaweb = get_metaweb()

SCENARIO_LABEL = Dict(
    "SSP126" => "SSP 1-2.6", 
    "SSP245" => "SSP 2-4.5",
    "SSP370" => "SSP 3-7.0"
)


# Find all species pairs that don't cooccur 
species = sort(collect(keys(sdms)))

bee_species = filter(contains("Bombus"), species)
plant_species = filter(!contains("Bombus"), species)

total_size = prod(size(sdms[species[begin]][:baseline][:range]))
overlap_mat = zeros(length(bee_species), length(plant_species))


for (i,b) in enumerate(bee_species), (j,p) in enumerate(plant_species)
    range_i = sdms[b][:baseline][:range]
    range_j = sdms[p][:baseline][:range]
    overlap_mat[i,j] = sum(range_i .& range_j) / total_size
end


function get_new_cooccurrence(sdms, overlap_mat, threshold, ssp, year)
    idx = findall(x -> x < threshold, overlap_mat)

    new_cooccurrence = Int64.(similar(first(sdms)[2][:baseline][:range]))
    new_cooccurrence.grid .= 0
    for i in idx
        b, p = bee_species[i[1]], plant_species[i[2]]

        range_i = sdms[b][:future][ssp][year][:range]
        range_j = sdms[p][:future][ssp][year][:range]
        
        new_overlap = range_i .& range_j
        new_cooccurrence += new_overlap
    end
    return new_cooccurrence
end



years = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]
ssps = ["SSP126", "SSP245", "SSP370"]

overlap_threshold = 0.05

gain_counts_by_ssp = Dict(
    [
        ssp=>[get_new_cooccurrence(sdms, overlap_mat, overlap_threshold, ssp, yr) for yr in years]
        for ssp in ssps
    ]
) 

extrema.(gain_counts_by_ssp[ssps[3]])


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

N_QUANTILES = 5


ref = gain_counts_by_ssp["SSP245"][end]
q_gain, qs = quantize_with_reference(gain_counts_by_ssp["SSP245"], ref; nbreaks=N_QUANTILES)

qs = quantile(values(nodata(ref, 0)), [i for i in LinRange(0,1,nbreaks+1)[2:end-1]])

begin 
f = Figure(size=(1500, 1000))
g = GridLayout(f[1,1])
for (j, ssp) in enumerate(ssps) 
    q_gain, qs = quantize_with_reference(gain_counts_by_ssp[ssp], ref)
    for (i, yr) in enumerate(years)


        ax = Axis(
            g[j,i], 
            title= j == 1 ? yr : "", 
            ylabel = i == 1 ? SCENARIO_LABEL[ssp] : "",
            ylabelrotation=0,
            ylabelfont=:bold,
            ylabelsize=24,
            aspect=DataAspect(),
            titlealign=:left,
            titlesize=24,
        )
        heatmap!(ax, first(sdms)[2][:baseline][:range], colormap=[:grey70])
        heatmap!(ax, q_gain[i], colormap=:PuBu_5)
        lines!(ax, county_poly, color=:grey20, linewidth=0.25)
        lines!(ax, state_poly, color=:grey20, linewidth=1.)
    end
end

    Colorbar(
        g[:,5],
        colormap=cgrad(:PuBu_5, N_QUANTILES, categorical = true),
        ticks = (0.5:4.5, string.(Int.(vcat(1,qs...)))),
        colorrange=(0.5,5.5),
        width= 30,
        height= 400,
        label = "Total New Cooccurrences",
        flipaxis=false,
        labelsize = 16
    )
    f
end 

save(joinpath("plots", "novel_cooccurrence.png"), f)

