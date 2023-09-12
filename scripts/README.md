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

The later steps are designed to be run (in large part) on a SLURM based cluster.
Repredocuing the entire analysis on a local machine with limited cores/RAM may
prove difficult. 

1. Clean data (`01_clean_data.jl`)
2. Engineer features (`02_engineer_features`)
    - `01_vgae`
    - `02_feature_learning`
        - reps
            - Each directory also contains a `vis.jl`
      file for producing the supplemental visualizatoin figures for each
3. Interaction prediction (`03_interaction_prediction`)
    - This directory contains a folder which fits each classification model on
      each set of species representations. Each is designed to run on SLURM
      clusters as parallel job-arrays.  
        - `logistic`
        - `rf`
        - `brt`
        - `xgboost`
        - `ensemble`
4. Fit SDMs (`04_make_sdms`)
    - `make_sdms.jl`
5. Compute overlap (`05_compute_overlap`)
    - foo


