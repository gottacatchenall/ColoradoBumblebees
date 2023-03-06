
load_data() = BeeData(create_interaction_data()..., load_occurrence_data(), load_environmental_data(), load_cooccurence_data())
 

function load_occurrence_data()
    bee_occ_df = CSV.read(datadir("public", "occurrence", "bees.csv"), DataFrame)
    plant_occ_df = CSV.read(datadir("public", "occurrence", "plants.csv"), DataFrame)
    cols = [:species, :latitude, :longitude]
    vcat(bee_occ_df[!, cols], plant_occ_df[!, cols])
end

load_environmental_data() = CSV.read(datadir("public", "environment", "covariates.csv"), DataFrame)

load_cooccurence_data() = CSV.read(datadir("embargo", "cooccurence", "clean", "cooccurence.csv"), DataFrame)

