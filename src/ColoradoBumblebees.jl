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
    using GeoInterface

    abstract type Site end 
    struct PikesPeak <: Site end
    struct Gothic <: Site end
    struct ElkMeadows <: Site end
    struct Metaweb <: Site end
    export Site, PikesPeak, Gothic, ElkMeadows, Metaweb


    abstract type FeatureType end 
    struct Phylogenetic <: FeatureType end
    struct Environment <: FeatureType end
    struct Spatial <: FeatureType end
    struct Temporal <: FeatureType end
    struct RelativeAbundance <: FeatureType end 
    struct Structural <: FeatureType end 
    export FeatureType, Phylogenetic, Environment, Spatial, Temporal, RelativeAbundance, Structural
    
    
    include(srcdir(joinpath("cleaning", "clean_raw_interactions.jl")))
    export clean_interactions

    include(srcdir(joinpath("cleaning", "cleaned_interactions_to_networks.jl")))
    export create_interaction_data

end 