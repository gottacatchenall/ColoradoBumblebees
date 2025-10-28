using DataFrames
using CSV
using Statistics
using Turing
using Random
using Distributions
using JSON 
using MCMCChains

include(joinpath("..", "src", "io.jl"))
include(joinpath("..", "src", "networks.jl"))
include(joinpath("..", "src", "phenology.jl"))

function get_phenology(
    data_dir, 
    species; 
    startdate = Date(2025, 5, 1),
    enddate = Date(2025, 10, 1)
)
    gbif_df =  get_gbif_data(data_dir)
    taxa_df = get_taxa_df(data_dir)

    idx = findfirst(x->x.species_name==species, eachrow(taxa_df))
    isnothing(idx) && return nothing
    key = taxa_df.species_key[idx]
    filter!(x->x.speciesKey == key, gbif_df)
    
    dts = [Date(DateTime(replace(d, "Z" => ""))) for d in gbif_df.eventDate]
    doys = [(d - Date(Year(d).value, 1,1)).value + 1 for d in dts]

    startdoy = (startdate - Date(2025, 1, 1)).value 
    enddoy = (enddate - Date(2025, 1, 1)).value

    enddoy - startdoy + 1

    phen = zeros(enddoy - startdoy + 1)
    doy = collect(startdoy:enddoy)
    for i in doys
        arr_idx = i - startdoy + 1
        if arr_idx > 0 && arr_idx <= length(phen)
            phen[arr_idx] += 1
        end
    end
    doy, phen
end

function fit_phenology(
    data_dir,
    artifact_dir, 
    species; 
    max_k=3, 
    startdate = Date(2025, 5, 1),
    enddate = Date(2025, 10, 1),
    num_samples = 10_000,
    burn_in = 5_000
)
    x,y = get_phenology(data_dir, species, startdate=startdate, enddate=enddate)
    result = fit_gmm(
        x, 
        y, 
        max_k;
        num_samples = num_samples,
        burn_in = burn_in
    )
    mkpath(joinpath(artifact_dir, species))
    write_gmm(result, joinpath(artifact_dir, species, "phenology.json"))
end

artifact_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/artifacts" : "./artifacts"
data_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/data" : "./data"

job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

species = sort(get_species_list(data_dir))
fit_phenology(data_dir, artifact_dir, species[job_id]; max_k = 2)



"""
x,y = get_phenology(data_dir, "Linaria dalmatica")

scatter(x,y)
result = fit_gmm(
    x, 
    y, 
    3;
    num_samples = 10_000,
    burn_in = 5_000
)


result

f = Figure()
plot_gmm(f, (1,1), result; title="")
f
"""