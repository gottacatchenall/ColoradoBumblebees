The scripts folder contains all of the files you need to run to recreate the
analysis.

It is divided into two parts: the `prerequisites` directory, which perpares a
few data products required to run the rest of the analysis. The second part is
the set of `julia` scripts used to complete the analysis. Each of the scripts is
designed in order to be as short and straightforward as possible, with much of
the functionality and methods being used by each script contained in the `src`
directory.  

## The prerequisites (`01_prerequisites`)

- Preparing the Phylogeny using BEAST
- Preparing the CHELSA layers 


## Producing analysis results  (`02_analysis`)

The later steps are designed to be partially run on a SLURM based cluster. 

1. Clean data (`01_clean_data.jl`)
2. Engineer features (`02_engineer_features`)
    - `01_vgae`
    - `02_feature_learning`
3. Interaction prediction (`03_interaction_prediction`)
    - `logistic`
    - `rf`
    - `brt`
    - `xgboost`
4. Fit SDMs (`04_make_sdms`)
    - `make_sdms.jl`
5. Compute overlap (`05_compute_overlap`)
    - foo
6. Compute interaction richness (`06_interaction_richness`)
    - foo

## Visualization (`03_visualization`)

