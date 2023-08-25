#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --job-name=vgae-fit 
#SBATCH --output=slurm-vgae_fit-%A.%a.out 
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G      
#SBATCH --time=01:00:00         
#SBATCH --array=1-40


module load julia/1.8.5

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia classification_fit.jl
