module ColoradoBumblebees
    using Reexport

    @reexport using AbstractTrees
    @reexport using BSON
    @reexport using Clustering 
    @reexport using CSV
    @reexport using CUDA
    @reexport using Dates
    @reexport using DataFrames
    @reexport using DelimitedFiles
    @reexport using Distributions
    @reexport using DrWatson    
    @reexport using EcologicalNetworks
    @reexport using Flux
    @reexport using Graphs
    @reexport using JSON
    @reexport using LinearAlgebra
    @reexport using Mangal
    @reexport using MLJ
    @reexport using MultivariateStats
    @reexport using NewickTree
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
end 