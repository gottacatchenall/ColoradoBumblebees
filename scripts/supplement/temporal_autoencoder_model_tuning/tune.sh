#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G               # memory per node
#SBATCH --time=30:00

module load cuda
module load julia
module load cudnn 


export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/FluxGPU"

julia ae_tuning.jl