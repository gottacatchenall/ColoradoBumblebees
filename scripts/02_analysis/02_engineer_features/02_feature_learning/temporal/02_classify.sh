#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --job-name=temporal_batch_fits 
#SBATCH --output=slurm-temporal_batch_fits-%A.%a.out 
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G      
#SBATCH --time=01:00:00         
#SBATCH --array=1-120 


module load julia/1.8.5

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia fit_classifier.jl
