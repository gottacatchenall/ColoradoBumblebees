# 🐝🪻🕸️ Mapping the Rewiring of Colorado Bumble Bee Pollination Networks

This repository contains code associated with the paper _Projecting the spatial rewiring of the bumble bee pollination network of the Southern Rocky Mountains_.

The GBIF occurrence data is already cleaned and available in `data/gbif_clean.csv`, however GBIF can be queried again via `scripts/gbif_query.jl`. It requires a `.env` file with GBIF login credentials. 

The code for running the analysis and producing visualizations is written in the Julia language. Instructions for installing Julia can be found [here](https://julialang.org/downloads/). Specifically, Julia v1.11 is what was used for development.

## Fitting SDMs and projecting ranges

The SDMs and projections are produced via `fit_sdms.sh` run on a SLURM-based cluster.

Running these yourself requires adjusting a few variables to point to the correct location of the path to data. Specifically the variables

```bash
# path where outputs will be written
export ARTIFACT_DIR="/scratch/mcatchen/ColoradoBees/artifacts" 

# path to the data dir within this repo
export DATA_DIR="/scratch/mcatchen/ColoradoBees/data" 

# path to where WorldClim layers are stored, specifically with their paths ending in "_bio_XX.tif", where XX is the number of the layer (this is how they are provided by default).
export WORLDCLIM_DIR="/project/def-tpoisot/mcatchen/WorldClim" 
```

## Creating visualizations

The scripts used to create the visualizations in the paper and supplement are in the `viz` directory. 

Each script is associated with a single figure. The final figure in the manuscript (Figure 5; the Sankey diagram) is contained within it's own directory within `viz` that has its own `Project.toml` file, which is necessary to avoid dependency conflicts between the package used to make the Sankey visualization (SankeyMakie) and the packages used in the main repository.



