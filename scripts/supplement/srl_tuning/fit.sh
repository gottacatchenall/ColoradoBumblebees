#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --mem=16G               # memory per node
#SBATCH --time=30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1


module load julia/1.8.5

export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/COBees"
export CLUSTER="true"

julia test_distributed_classifier_fit.jl
