#!/bin/bash
#SBATCH --job-name=temporal_batch_fits 
#SBATCH --output=slurm-%A.%a.out 
#SBATCH --nodes=1               
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=16G      
#SBATCH --time=00:20:00         
#SBATCH --array=1-64 


module load julia/1.8.5

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia test_distributed_classifier_fit.jl
