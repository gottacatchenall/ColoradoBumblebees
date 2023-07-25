function load_data()
    return BeeData(
        load_interaction_data()...,
        load_occurrence_data(),
        load_environmental_data(),
    )
end

load_occurrence_data() = CSV.read(datadir("public", "occurrence", "occurrence.csv"), DataFrame)

function load_environmental_data()
    return CSV.read(datadir("public", "environment", "covariates.csv"), DataFrame)
end


