
include(joinpath("..", "src", "sdm.jl"))
include(joinpath("..", "src", "networks.jl"))

function main(
    data_dir, 
    artifact_dir, 
    species_name;
    cluster = false,
    k = 5,
    class_balance = 1.
)

    bbox = (left=-109.7, right=-101.8, bottom=34.5,top=42.5)

    #cluster ? 
    #L = [SDMLayer(RasterData(CHELSA2, BioClim); layer=i, bbox...) for i in 1:19]
    #L = [Float32.(l) for l in L]

    cluster_chelsa_dir = "/home/mcatchen/projects/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees/SimpleSDMDatasets/CHELSA2/BioClim/"
    chelsa_paths = [joinpath(cluster_chelsa_dir, findfirst(isequal(x), readdir(cluster_chelsa_dir))) for x in ["_bio$(i)_" for i in 1:19]]
    L = [SDMLayer(p; bbox...) for p in chelsa_paths]
    L = [Float32.(l) for l in L]

    occ_records = get_occurrences(data_dir)
    occs = split_occurrences_into_species(occ_records)
    occ = occs[species_name]

    pred, unc, fit_stats, pr, ab = fit_sdm(
        occ, 
        L; 
        k = k,
        class_balance = class_balance
    )

    write_sdm_artifacts(
        artifact_dir, 
        species_name, 
        Dict(
            :prediction => pred,
            :uncertainty => unc,
            :presences => pr,
            :absences => ab, 
            :metrics => fit_stats
        )
    )
end 

artifact_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/artifacts" : "./artifacts"
data_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/data" : "./data"

job_id = ENV["SLURM_ARRAY_TASK_ID"]

species = sort(get_species_list(data_dir))
main(data_dir, artifact_dir, species[job_id], cluster = "cluster" in ARGS)

