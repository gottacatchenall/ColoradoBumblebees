#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=SDMs 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-180 
#SBATCH --time=12:00:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

module load julia/1.11.3
srun --unbuffered julia -e '
    include(joinpath("..", "src", "sdms.jl"))
    include(joinpath("..", "src", "networks.jl"))

    artifact_dir = "/scratch/mcatchen/ColoradoBees/artifacts" 
    data_dir = "/scratch/mcatchen/ColoradoBees/data"
    worldclim_dir = "/project/def-tpoisot/mcatchen/WorldClim" 

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species_list = sort(get_species_list(data_dir))
    species_name = species_list[job_id]

    hyperparams = get_hyperparameters(data_dir, species_name)
    create_species_distribution_models(data_dir, artifact_dir, worldclim_dir, species_name; hyperparams...)
'