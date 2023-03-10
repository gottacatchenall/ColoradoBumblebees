using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

using SpeciesDistributionToolkit

_scenario_to_projection(scenario, model=GFDL_ESM4) = Projection(scenario.ssp, model)
function _scenario_to_year_pair(scenario)
    return scenario.years.startyear => scenario.years.endyear
end


for s in scenarios()[2:end]
    tmp = [
        SimpleSDMPredictor(
            RasterData(CHELSA2, BioClim),
            _scenario_to_projection(s);
            timespan=_scenario_to_year_pair(s),
            layer=l,
            EXTENT...,
        ) for l in BIOLAYERS
    ]
end