#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=TuningSDMs 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-180 
#SBATCH --time=12:00:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

export ARTIFACT_DIR="/scratch/mcatchen/ColoradoBees/artifacts"
export DATA_DIR="/scratch/mcatchen/ColoradoBees/data"
export WORLDCLIM_DIR="/project/def-tpoisot/mcatchen/WorldClim"

module load julia/1.11.3
srun --unbuffered julia -e '
    include(joinpath("..", "src", "sdms.jl"))
    include(joinpath("..", "src", "networks.jl"))

    artifact_dir = ENV["ARTIFACT_DIR"]
    data_dir = ENV["DATA_DIR"]
    worldclim_dir = ENV["WORLDCLIM_DIR"]

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species = sort(get_species_list(data_dir))
    
    tune_hyperparameters(
        data_dir, 
        artifact_dir, 
        worldclim_dir,
        species[job_id];
        k = 5,
        class_balances = 0.5:0.5:3,
        pseudoabsence_buffer_distances = 5.0:5.0:25,
        max_depths = 4:2:10,
    )
'