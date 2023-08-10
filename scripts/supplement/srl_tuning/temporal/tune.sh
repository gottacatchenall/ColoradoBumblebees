#!/bin/bash
#SBATCH --account=def-gonzalez
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G               # memory per node
#SBATCH --time=1:06:00

module load cuda
module load julia
module load cudnn 


export JULIA_DEPOT_PATH="/project/def-gonzalez/mcatchen/JuliaEnvironments/FluxGPU"

julia ae_tuning.jl