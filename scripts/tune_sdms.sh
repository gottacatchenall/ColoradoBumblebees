#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=TuningSDMs 
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
    include(joinpath("..", "src", "sdms.jl"))
    include(joinpath("..", "src", "networks.jl"))

    artifact_dir = "/scratch/mcatchen/ColoradoBees/artifacts" 
    data_dir = "/scratch/mcatchen/ColoradoBees/data"
    worldclim_dir = "/project/def-tpoisot/mcatchen/WorldClim" 

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species = sort(get_species_list(data_dir))
    
    tune_hyperparameters(
        data_dir, 
        artifact_dir, 
        worldclim_dir,
        species_name;
        k = 5,
        class_balances = 0.5:0.5:3,
        pseudoabsence_buffer_distances = 5.0:5.0:25,
        max_depths = 4:2:10,
    )
'