
include(srcdir("exports.jl"))
include(srcdir("constants.jl"))

include(srcdir("types", "sites.jl"))
include(srcdir("types", "species.jl"))
include(srcdir("types", "data.jl"))
include(srcdir("types", "embedding.jl"))
include(srcdir("types", "representations.jl"))
include(srcdir("types", "classification_models.jl"))
include(srcdir("types", "classification_fits.jl"))
include(srcdir("types", "confusion_matrix.jl"))
include(srcdir("types", "species_distributions.jl"))

include(srcdir("load_data", "load_chelsa.jl"))
include(srcdir("load_data", "load_networks.jl"))
include(srcdir("load_data", "load_interactions.jl"))
include(srcdir("load_data", "load_data.jl"))
include(srcdir("load_data", "load_phylogeny.jl"))
include(srcdir("load_data", "load_phenology.jl"))

include(srcdir("data_cleaning", "clean_interactions.jl"))
include(srcdir("data_cleaning", "clean_environment_covariates.jl"))


# Species Representation Learning 
include(srcdir("species_representation_learning", "structural", "shared.jl"))
include(srcdir("species_representation_learning", "structural", "svd.jl"))
include(srcdir("species_representation_learning", "structural", "lfsvd.jl"))

include(srcdir("species_representation_learning", "node2vec.jl"))

include(srcdir("species_representation_learning", "phylogenetic", "node2vec.jl"))
include(srcdir("species_representation_learning", "phylogenetic", "simulated_traits.jl"))

include(srcdir("species_representation_learning", "environment", "kmeans.jl"))

include(srcdir("species_representation_learning", "relative_abundance.jl"))

include(srcdir("species_representation_learning", "temporal", "dense.jl"))
include(srcdir("species_representation_learning", "temporal", "recurrent.jl"))

include(srcdir("species_representation_learning", "feature_dataframe.jl"))
include(srcdir("species_representation_learning", "embed.jl"))

# Classification model fit
include(srcdir("classification_models", "classification_stats.jl"))
include(srcdir("classification_models", "balance_sample.jl"))
include(srcdir("classification_models", "fit_classifier.jl"))
include(srcdir("classification_models", "batch_fit.jl"))

include(srcdir("artifacts", "shared.jl"))
include(srcdir("artifacts", "representation_artifacts.jl"))
include(srcdir("artifacts", "classification_artifacts.jl"))

# SDMs

include(srcdir("species_distributions", "distribution_models.jl"))
include(srcdir("species_distributions", "utils.jl"))
include(srcdir("species_distributions", "fit_distribution.jl"))
