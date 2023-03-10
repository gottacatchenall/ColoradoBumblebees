using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

using SpeciesDistributionToolkit

for s in scenarios()
    tmp = [
        SimpleSDMPredictor(
            RasterData(CHELSA2, BioClim),
            _scenario_to_projection(scenario);
            timespan=_scenario_to_year_pair(scenario),
            layer=l,
            EXTENT...,
        ) for l in BIOLAYERS
    ]
end