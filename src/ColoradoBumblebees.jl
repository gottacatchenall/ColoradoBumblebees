#= 

Structure

- Clean and prepare data
    - 1. Convert raw interaction CSVs to cleaned CSVs
    - 2. Convert cleaned CSVs into Mangal data structures
    
- Build features: create a dataframe of each species' features and each way
   they can be constructed 

    Features:

        1. Temporal
        --------------------------------------------------

            Data: 
                - 1. Interaction-derived temporal distributions across years/plots
                - 2. (Maybe) phenology timeseries from RMBL

            Methods for feature engineering:
                - Fit a Gaussian and report mean and variance
                    - Do this stratified across spatial sites or pooled
                - (RNN?)  Autoencoder trained on gaussian phenology abundance data (or the real data)

        2. Phylogenetic
        --------------------------------------------------

            Data:
                - Estimated phylogeny

            Methods for feature engineering:
                - Simulating N traits across a rate spectrum and PCAing
            
            
        3. **Spatial** 
        --------------------------------------------------

            Data: 
                - binary occurrence rasters over the extent

            Methods for feature engineering:
                - autoencoders

        4. **Environmental**
        --------------------------------------------------

            Data:
                - Environment vectors for each occurrence location

            Methods for feature engineering:
                - Fit MV Normal
                - Fit K-means to output of center locations
                - Fit autoencoder on batches of environmental vectors
                

        5. **Relative abundance** 
        --------------------------------------------------
            Data: 
                - RA from interaction data
                - (Possibly): RA from phenology data

            Methods for feature engineering:
                

        6. **Structural**
        --------------------------------------------------
            Data:
                - Interaction matrix / Occurrence matrix

            Methods for feature engineering:
                - PCA/SVD
=#

# end

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
    using Distributions
    #using MLJ
    using Clustering
    using Random
    using ProgressMeter
    using ParameterSchedulers
    using ParameterSchedulers: Scheduler
    
    const extent = (bottom=34., top=44., left=-110.5, right=-103.5)
    export extent 


    include(srcdir("util.jl"))

    export bee, plant, pollinator, bee
    export plants, bees, pollinators

    
    abstract type Site end 
    struct PikesPeak <: Site end
    struct Gothic <: Site end
    struct ElkMeadows <: Site end
    struct Metaweb <: Site end
    export Site, PikesPeak, Gothic, ElkMeadows, Metaweb, sitename
    sitename(::Type{PikesPeak}) = "Pikes Peak"
    sitename(::Type{Gothic}) = "Gothic"
    sitename(::Type{ElkMeadows}) = "Elk Meadows"


    struct Bee 
        name
        mangalnode::MangalNode
    end

    struct Plant
        name
        mangalnode::MangalNode
    end 

    Base.show(io::IO, plant::Plant) = Base.show(io, "🌷 $(plant.name)")
    Base.show(io::IO, bee::Bee) = Base.show(io, "🐝 $(bee.name)")
    export Bee, Plant


    struct Interaction{T<:Site}
        bee::Bee
        plant::Plant
        int::MangalInteraction
        elevation 
        time
    end
    Base.show(io::IO, int::Interaction{T}) where T = Base.print(io,"🐝 $(int.bee.name) ↔️ 🌷 $(int.plant.name) at $(sitename(T)) on $(monthname(int.time)) $(day(int.time)), $(year(int.time))")
    

    struct BeeData
        bees::Vector{Bee}
        plants::Vector{Plant}
        interactions::Vector{Interaction}
        occurrence::DataFrame
        environment::DataFrame
    end 
    Base.show(io::IO, bd::BeeData) = Base.print(io, "Pollination dataset with $(length(bd.interactions)) interactions")
 
    interactions(bd::BeeData) = bd.interactions

    function interactions(bd::BeeData, sp1, sp2)
        thisbee, thisplant = split(sp1.name, " ")[1] == "Bombus" ? (sp1, sp2) : (sp2, sp1)
        filter(int ->  bee(int) == thisbee && plant(int) == thisplant, interactions(bd))
    end
    
    function interactions(bd::BeeData, sp1)
        isbee = split(sp1.name, " ")[1] == "Bombus"

        otherspecies = isbee ? plants(bd) : bees(bd)
        Is = []
        for sp2 in otherspecies
            thisbee = isbee ? sp1 : sp2
            thisplant = isbee ? sp2 : sp1
            Is = vcat(Is..., findall(int ->  bee(int) == thisbee && plant(int) == thisplant, interactions(bd))...)
        end
        interactions(bd)[Is]
    end


    occurrence(bd::BeeData) = bd.occurrence
    environment(bd::BeeData) = bd.environment
    
    plant(int::MangalInteraction) = int.to
    plant(int::Interaction) = int.plant
    plant(bd::BeeData, str) = plants(bd)[findfirst(x->x.name==str, plants(bd))]


    bee(int::MangalInteraction) = int.from
    bee(int::Interaction) = int.bee
    bee(bd::BeeData, str) = bees(bd)[findfirst(x->x.name==str,bees(bd))]

    bees(data) = data.bees
    plants(data) = data.plants
    


    function phenology(data)    
        species = vcat(bees(data)..., plants(data)...)
        firstdoy, lastdoy = extrema([dayofyear(i.time) for i in interactions(data)])

        dict = Dict()
        for sp in species
            abundances = zeros(Int32, lastdoy-firstdoy+1)
            ints = interactions(data, sp)
            for i in ints
                abundances[dayofyear(i.time) - firstdoy + 1] += 1
            end
            merge!(dict, Dict(sp=>abundances))
        end
        dict
    end

    export BeeData, Interaction, interactions, occurrence, environment, phenology

    function balance_sample(y, I, batch_size=64, true_pct=0.5)
        itrue = findall(x-> x == true, y)
        ifalse = findall(x-> x == false, y)
    
        filter!(x->x ∈ I, itrue)
        filter!(x->x ∈ I, ifalse)
    
        ntrue, nfalse = Int32(floor(batch_size*true_pct)),Int32(floor(batch_size*(1-true_pct)))
        [shuffle(itrue)[1:ntrue]..., shuffle(ifalse)[1:nfalse]...]
    end
    export balance_sample


    abstract type FeatureType end 
    abstract type Phylogenetic <: FeatureType end
    abstract type Environment <: FeatureType end
    abstract type Spatial <: FeatureType end
    abstract type Temporal <: FeatureType end
    abstract type RelativeAbundance <: FeatureType end 
    abstract type Structural <: FeatureType end 
    export FeatureType, Phylogenetic, Environment, Spatial, Temporal, RelativeAbundance, Structural
    
    export getfeatures
    export load_data, load_occurrence_data

    include(srcdir(joinpath("cleaning", "clean_raw_interactions.jl")))
    export clean_interactions

    include(srcdir(joinpath("cleaning", "create_networks.jl")))
    export create_interaction_data

    include(srcdir(joinpath("cleaning", "clean_environmental_covariates.jl")))
    export create_environmental_covariate_data

    include(srcdir(joinpath("cleaning", "load_data.jl")))
    export load_data


    export outdim
    include(srcdir("features", "spatial.jl"))
    export KMeansSpatialEmbedding

    include(srcdir("features", "environment", "kmeans.jl"))
    include(srcdir("features", "environment", "autoencoder.jl"))

    export KMeansEnvironmentEmbedding, EnvironmentAutoencoder


    abstract type AutoencoderType end 
    struct Standard <: AutoencoderType end
    struct Variational <: AutoencoderType end 

    include(srcdir("features", "temporal.jl"))
    export TemporalAutoencoder, Variational, Standard

    const TEMPORAL_INPUT_DIM = 147
    export TEMPORAL_INPUT_DIM


    include(srcdir("models", "sandbox.jl"))
    export feature_dataframe, label_dataframe

    include(srcdir("models", "confusionmatrix.jl"))
    export computemeasures

end 