#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=FinishTuningSDMs 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=[1,129,19,73,103,178,46,38,92,110,29,135,93,172,14,151] 
#SBATCH --time=8:00:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

module load julia/1.11.3
julia -e '
    include(joinpath("..", "src", "sdms.jl"))
    include(joinpath("..", "src", "networks.jl"))

    artifact_dir = "/scratch/mcatchen/ColoradoBees/artifacts" 
    data_dir = "/scratch/mcatchen/ColoradoBees/data"
    worldclim_dir = "/project/def-tpoisot/mcatchen/WorldClim" 

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species = sort(get_species_list(data_dir))

    finish_tuning_hyperparameters(
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