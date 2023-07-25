module ColoradoBumblebees
    using AbstractTrees
    using CSV
    using Dates
    using DataFrames
    using DelimitedFiles
    using Distributions
    using DrWatson    
    using EcologicalNetworks
    using Flux
    using Graphs
    using JSON
    using LinearAlgebra
    using Mangal
    using MLJ
    using MultivariateStats
    using NewickTree
    using Phylo
    using ProgressMeter
    using Random
    using SpeciesDistributionToolkit
    using Statistics
    using StatsBase
    using Word2Vec
        
    
    # Constants
    const EXTENT = (bottom=34.0, top=44.0, left=-110.5, right=-103.5)
    const TEMPORAL_INPUT_DIM = 147
    const BIOLAYERS = ["BIO$i" for i in 1:19]

    export TEMPORAL_INPUT_DIM, EXTENT, BIOLAYERS


    # -----------------------------------------------------------------
    # Types 
     
    # Sites
    export Site
    export PikesPeak, Gothic, ElkMeadows

    # Species
    export Bee, Plant

    # Data
    export BeeData
    export Interaction

    # Features
    export FeatureType

    export Phylogenetic
    export Environment
    export Spatial
    export Temporal
    export RelativeAbundance
    export Structural 


    # -----------------------------------------------------------------
    # Methods 

    # Getters and utils
    export plant
    export bee
    export plants
    export bees

    export cooccurence
    export occurrence
    export interactions
    export phenology
    export metaweb
    export environment
    export sitename

    # Load data
    export load_data 
    export load_occurrence_data
    export load_baseline_layers
    export load_interaction_data
    export load_newick
    export load_phenology

    # Cleaning 
    export clean_interactions
    export clean_environmental_covariate_data
    
    # Features
    export outdims
    export getfeatures 
    export feature_dataframe, label_dataframe

    export MetawebSVD, LFSVD 
    export SimulatedTraits, PhylogeneticNode2Vec
    export KMeansEnvironmentEmbedding
    export Pooled
    export DenseAutoencoder, Standard, Variational
    export RecurrentAutoencoder
    
    # Feature Learning

    export crossvalidation
    export batch
    export fit_model
    export computemeasures

    # -----------------------------------------------------------------
    # Includes 


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

    include(srcdir("feature_learning", "structural", "svd.jl"))
    include(srcdir("feature_learning", "structural", "lfsvd.jl"))

    include(srcdir("feature_learning", "node2vec.jl"))
    include(srcdir("feature_learning", "phylogenetic", "node2vec.jl"))
    include(srcdir("feature_learning", "phylogenetic", "simulated_traits.jl"))

    include(srcdir("feature_learning", "environment", "kmeans.jl"))

    include(srcdir("feature_learning", "relative_abundance.jl"))

    include(srcdir("feature_learning", "temporal", "dense.jl"))
    include(srcdir("feature_learning", "temporal", "recurrent.jl"))

    include(srcdir("feature_learning", "feature_dataframe.jl"))
    include(srcdir("feature_learning", "confusion_matrix.jl"))
    include(srcdir("feature_learning", "balance_sample.jl"))
    include(srcdir("feature_learning", "crossvalidation.jl"))
    include(srcdir("feature_learning", "batch_fit.jl"))

end 