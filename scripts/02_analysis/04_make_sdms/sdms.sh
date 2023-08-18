#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --job-name=make_sdms 
#SBATCH --output=slurm-make_sdms.%A.%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-180 
#SBATCH --time=45:00         


module load julia/1.8.5

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia make_sdms.jl
