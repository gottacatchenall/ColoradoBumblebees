
function load_data()
    return BeeData(
        create_interaction_data()...,
        load_occurrence_data(),
        load_environmental_data(),
        load_cooccurence_data(),
    )
end

function load_occurrence_data()
    bee_occ_df = CSV.read(datadir("public", "occurrence", "bees.csv"), DataFrame)
    plant_occ_df = CSV.read(datadir("public", "occurrence", "plants.csv"), DataFrame)
    cols = [:species, :latitude, :longitude]
    return vcat(bee_occ_df[!, cols], plant_occ_df[!, cols])
end

function load_environmental_data()
    return CSV.read(datadir("public", "environment", "covariates.csv"), DataFrame)
end

function load_cooccurence_data()
    return CSV.read(
        datadir("embargo", "cooccurence", "clean", "cooccurence.csv"), DataFrame
    )
end
