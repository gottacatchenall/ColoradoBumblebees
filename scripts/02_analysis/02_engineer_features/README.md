`02_engineer_features`

First run the scripts in `01_vgae`. This runs the graph-autoencoder (GAE) models
which are implemented in `python` using `PyTorch` and `pytorch-geometric`. These
scripts produce files in `artifacts/vgae` which are used in the next step. 

Then run the scripts in `02_feature_learning`, which generate
species-representations using each model across a variety of hyperparameter
configurations. 

There are also scripts for visualizing the difference between
hyperparameterizations of the same model, which are included in the manuscript
supplement. 
