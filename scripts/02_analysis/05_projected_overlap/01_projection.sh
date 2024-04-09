#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --job-name=project_overlap 
#SBATCH --output=slurm-project_overlap.%A.%a.out
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-25
#SBATCH --time=45:00         



export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

module load julia/1.10.0

julia project_future_overlap.jl