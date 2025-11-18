#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=Phenologies 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-180 
#SBATCH --time=3:00:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

module load julia/1.11.3
julia -e '    
    include(joinpath("..", "src", "io.jl"))
    include(joinpath("..", "src", "networks.jl"))
    include(joinpath("..", "src", "phenology.jl"))

    artifact_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/artifacts" : "./artifacts"
    data_dir = "cluster" in ARGS ? "/scratch/mcatchen/ColoradoBees/data" : "./data"

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species = sort(get_species_list(data_dir))
    fit_phenology(data_dir, artifact_dir, species[job_id]; max_k = 3)
'