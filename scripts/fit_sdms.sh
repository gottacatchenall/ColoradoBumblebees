#!/bin/bash
#SBATCH --account=def-tpoisot
#SBATCH --job-name=SDMs 
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-197 
#SBATCH --time=0:30:00         

export JULIA_DEPOT_PATH="/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees"

module load julia/1.11.3
julia -e '
    include(joinpath("..", "src", "sdms.jl"))
    include(joinpath("..", "src", "networks.jl"))

    artifact_dir = "/scratch/mcatchen/ColoradoBees/artifacts" 
    data_dir = "/scratch/mcatchen/ColoradoBees/data"
    chelsa_dir = "/home/mcatchen/scratch/ColoradoBees/data/CHELSA" 

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    species = sort(get_species_list(data_dir))
    create_sdms(data_dir, artifact_dir, chelsa_dir, species[job_id])
'