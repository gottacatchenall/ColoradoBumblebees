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

function get_phenology(data_dir, species)
    gbif_df =  get_gbif_data(data_dir)
    taxa_df = get_taxa_df(data_dir)

    idx = findfirst(x->x.species_name==species, eachrow(taxa_df))
    isnothing(idx) && return nothing
    key = taxa_df.species_key[idx]
    filter!(x->x.speciesKey == key, gbif_df)
    
    dts = [Date(DateTime(replace(d, "Z" => ""))) for d in gbif_df.eventDate]
    
    doys = [(d - Date(Year(d).value, 1,1)).value + 1 for d in dts]

    phen = zeros(366)
    doy = collect(1:366)
    for i in doys
        phen[i] += 1
    end
    doy, phen
end

function fit_phenology(data_dir, artifact_dir, species; max_k=3)
    x,y = get_phenology(data_dir, species)
    result = fit_gmm(
        x, 
        y, 
        max_k;
    )
    mkpath(joinpath(artifact_dir, species))
    write_gmm(result, joinpath(artifact_dir, species, "phenology.json"))
end

artifact_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/artifacts" : "./artifacts"
data_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/data" : "./data"

job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

species = sort(get_species_list(data_dir))
fit_phenology(data_dir, artifact_dir, species[job_id])



for sp in species
    fit_phenology(data_dir, artifact_dir, sp)
end 

