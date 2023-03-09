module ColoradoBumblebees
    using DrWatson
    using DataFrames, CSV
    using Dates
    using Mangal
    using Statistics
    using GeoInterface
    using SpeciesDistributionToolkit
    using MultivariateStats
    using Flux
    using LinearAlgebra
    using Distributions
    using MLJ
    using Clustering
    using Random
    using ProgressMeter
    using ParameterSchedulers
    using ParameterSchedulers: Scheduler
    using EcologicalNetworks
    using NewickTree
    using AbstractTrees
    using Phylo

    # _______________________________________________________
    #
    # Exports
    #
    # _______________________________________________________
    
    # __________________
    # Constants
    # __________________
    const EXTENT = (bottom=34., top=44., left=-110.5, right=-103.5)
    const TEMPORAL_INPUT_DIM = 147
    export TEMPORAL_INPUT_DIM, EXTENT

    # __________________
    # Types
    # __________________

    export Site, PikesPeak, Gothic, ElkMeadows, Metaweb
    export Bee, Plant
    
    # Features
    export FeatureType, Phylogenetic, Environment, Spatial, Temporal, RelativeAbundance, Structural

    # Data
    export BeeData
    export Interaction

    # __________________
    # Methods
    # __________________
    
    # Getters and utils
    export cooccurence
    export occurrence
    export interactions
    export phenology
    export metaweb
    export environment
    export sitename

    export bee, pollinator
    export plant
    export bees, pollinators
    export plants

    # Data cleaning
    export clean_interactions
    export create_interaction_data
    export create_environmental_covariate_data
    export create_cooccurence_data

    # Models
    export feature_dataframe, label_dataframe
    export computemeasures, computemeasures_mlj
    export getfeatures
    export balance_sample

    # Data loading
    export load_data
    export load_occurrence_data
    export load_cooccurence_data
    export load_environmental_data

    # Features
    export KMeansSpatialEmbedding
    export KMeansEnvironmentEmbedding
    export EnvironmentAutoencoder
    export CooccurencePCA
    export MetawebSVD
    export Autoencoder
    export Variational
    export Standard
    export Pooled 
    export Hierarchical

    include(srcdir("types", "sites.jl"))
    include(srcdir("types", "species.jl"))
    include(srcdir("types", "data.jl"))
    include(srcdir("types", "features.jl"))

    include(srcdir("utils", "cooccurence.jl"))
    include(srcdir("utils", "interactions.jl"))
    include(srcdir("utils", "phenology.jl"))
    include(srcdir("utils", "load_data.jl"))

    include(srcdir("cleaning", "clean_raw_interactions.jl"))
    include(srcdir("cleaning", "create_networks.jl"))
    include(srcdir("cleaning", "clean_environmental_covariates.jl"))
    include(srcdir("cleaning", "create_cooccurence.jl"))

    include(srcdir("features", "spatial.jl"))
    include(srcdir("features", "environment", "kmeans.jl"))
    include(srcdir("features", "environment", "autoencoder.jl"))
    include(srcdir("features", "structural", "svd.jl"))
    include(srcdir("features", "temporal.jl"))
    include(srcdir("features", "relativeabundance.jl"))
    include(srcdir("features", "phylogenetic", "hierarchical.jl"))
    include(srcdir("features", "phylogenetic", "node2vec.jl"))
    include(srcdir("features", "phylogenetic", "simulated_traits.jl"))

    include(srcdir("models", "dataframes.jl"))
    include(srcdir("models", "confusionmatrix.jl"))
    include(srcdir("models", "balance_sample.jl"))

end 