sdmdir() = ColoradoBumblebees.CLUSTER ? joinpath("/scratch/mcatchen/BeeSDMs") : joinpath(artifactdir(), "sdms") 
sdmdir(sp::S) where S<:Species = joinpath(sdmdir(), sp.name)
sdmdir(sdm::SpeciesDistribution) = joinpath(sdmdir(sdm.species), string(sdm.timespan), string(sdm.scenario))

function ColoradoBumblebees.save(sdm::SpeciesDistribution)
    sdm_dir = sdmdir(sdm)
    run(`mkdir -p $sdm_dir`)

    SpeciesDistributionToolkit.save(joinpath(sdm_dir, "prediction.tif"), sdm.probability)
    SpeciesDistributionToolkit.save(joinpath(sdm_dir, "uncertainty.tif"), sdm.uncertainty)
    json_string = JSON.json(sdm.fit_stats)
    open(joinpath(sdm_dir, "fit.json"), "w") do f
        JSON.print(f, json_string)
    end
end

_string_to_scenario(x) = begin
    Dict(
        "Baseline" => Baseline,
        "SSP1_26" => SSP1_26,
        "SSP2_45" => SSP2_45,
        "SSP3_70" => SSP3_70
    )[x]
end

function load_sdm(sp::Species, timespan::Type{T}, scenario::Type{S}) where {T<:Timespan, S<:Scenario}
    _load_sdm(joinpath(sdmdir(sp), string(timespan), string(scenario)))
end 

function _load_sdm(path)
    data = load_data()
    pred = SimpleSDMPredictor(joinpath(path, "prediction.tif"))
    uncert = SimpleSDMPredictor(joinpath(path, "uncertainty.tif"))
    
    fit_stats = JSON.parse(JSON.parsefile(joinpath(path, "fit.json")))

    sp_string, yr_string, scenario_string = split(path, "/")[end-2:end]
    yrs = Year.(parse.(Int32,split(yr_string, "_")))

    species = contains(sp_string, "Bombus") ? bee(data, sp_string) : plant(data, sp_string)

    scenario = _string_to_scenario(scenario_string)
    timespan = ColoradoBumblebees.Timespan{yrs[1],yrs[2]}

    SpeciesDistribution(species, pred, uncert, fit_stats, timespan, scenario)
end

