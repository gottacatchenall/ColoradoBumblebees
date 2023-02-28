#!/bin/bash


#SBATCH --time=12:00:00
#SBATCH --account=def-gonzalez
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

export JULIA_DEPOT_PATH = "/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"

module load julia/1.8.1
echo "Launching script..."
julia -t 64 intro.jl
