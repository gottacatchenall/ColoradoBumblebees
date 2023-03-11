#!/bin/bash


#SBATCH --time=24:00:00
#SBATCH --account=def-gonzalez
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128G
#SBATCH --ntasks=203



export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"

module load julia/1.8.5
echo "Launching script..."
julia -t 203 make_sdms.jl
# julia -t 6 make_sdms.jl
