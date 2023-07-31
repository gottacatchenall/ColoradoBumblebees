
include(srcdir("exports.jl"))
include(srcdir("constants.jl"))

include(srcdir("types", "sites.jl"))
include(srcdir("types", "species.jl"))
include(srcdir("types", "data.jl"))
include(srcdir("types", "features.jl"))

include(srcdir("load_data", "load_chelsa.jl"))
include(srcdir("load_data", "load_networks.jl"))
include(srcdir("load_data", "load_interactions.jl"))
include(srcdir("load_data", "load_data.jl"))
include(srcdir("load_data", "load_phylogeny.jl"))
include(srcdir("load_data", "load_phenology.jl"))

include(srcdir("data_cleaning", "clean_interactions.jl"))
include(srcdir("data_cleaning", "clean_environment_covariates.jl"))


# Species Representation Learning 
include(srcdir("feature_learning", "structural", "shared.jl"))
include(srcdir("feature_learning", "structural", "svd.jl"))
include(srcdir("feature_learning", "structural", "lfsvd.jl"))

include(srcdir("feature_learning", "node2vec.jl"))

include(srcdir("feature_learning", "phylogenetic", "node2vec.jl"))
include(srcdir("feature_learning", "phylogenetic", "simulated_traits.jl"))

include(srcdir("feature_learning", "environment", "kmeans.jl"))

include(srcdir("feature_learning", "relative_abundance.jl"))

include(srcdir("feature_learning", "temporal", "dense.jl"))
include(srcdir("feature_learning", "temporal", "recurrent.jl"))

# Classification model fit
include(srcdir("classification_models", "feature_dataframe.jl"))
include(srcdir("classification_models", "confusion_matrix.jl"))
include(srcdir("classification_models", "balance_sample.jl"))
include(srcdir("classification_models", "crossvalidation.jl"))
include(srcdir("classification_models", "batch_fit.jl"))
