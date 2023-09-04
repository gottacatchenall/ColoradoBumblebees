module ColoradoBumblebees
    using Reexport

    @reexport using AbstractTrees
    @reexport using BSON
    @reexport using Clustering 
    @reexport using Combinatorics
    @reexport using CSV
    @reexport using CUDA
    @reexport using Dates
    @reexport using DataFrames
    @reexport using DelimitedFiles
    @reexport using Distributions
    @reexport using DrWatson    
    @reexport using EcologicalNetworks
    @reexport using EvoTrees
    @reexport using Flux
    @reexport using Graphs
    @reexport using JSON
    @reexport using LinearAlgebra
    @reexport using Mangal
    @reexport using MLJ
    @reexport using MultivariateStats
    @reexport using NewickTree
    @reexport using OffsetArrays
    @reexport using Phylo
    @reexport using ProgressMeter
    @reexport using Random
    @reexport using SpeciesDistributionToolkit
    @reexport using Statistics
    @reexport using StatsBase
    @reexport using Word2Vec


    Random.seed!(3141592653589793)

    GPU_AVAILABLE = CUDA.has_cuda_gpu()
    CLUSTER = "CLUSTER" ∈ collect(keys(ENV))

    include(srcdir("includes.jl"))

    const BEST_REPRESENTATIONS = [
        PhylogeneticNode2Vec(number_of_walks=250, walk_length=100, embedding_dim=16),
        Pooled(),
        RecurrentAutoencoder{Standard}(
            rnn_dims = [1, 8, 1],
            encoder_dims = [TEMPORAL_INPUT_DIM, 16],
            decoder_dims = [16, TEMPORAL_INPUT_DIM],
            opt = ADAM(0.005)
        ),
        LFSVD(truncation_dims=6, embed_dims=6),
        KMeansEnvironmentEmbedding(k=7),
    ]
end 