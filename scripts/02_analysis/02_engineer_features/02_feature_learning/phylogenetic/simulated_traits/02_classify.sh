#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --job-name=simulated_traits_fit 
#SBATCH --output=slurm-simulated_traits_fit-%A.%a.out 
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G      
#SBATCH --time=01:00:00         
#SBATCH --array=1-24

module load cuda
module load julia/1.8.5
module load cudnn 

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia kmeans_tuning.jl
