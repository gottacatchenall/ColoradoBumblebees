using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees
using SpeciesDistributionToolkit

_scenario_to_projection(scenario, model=GFDL_ESM4) = Projection(scenario.ssp, model)
function _scenario_to_year_pair(scenario)
    return scenario.years.startyear => scenario.years.endyear
end

N_ATTEMPTS = 15

for s in scenarios()[2:end]
    @info s
    for l in BIOLAYERS
        @info l
        ct = 0
        while ct < N_ATTEMPTS
            try
                tmp = 
                    SimpleSDMPredictor(
                        RasterData(CHELSA2, BioClim),
                        _scenario_to_projection(s);
                        timespan=_scenario_to_year_pair(s),
                        layer=l,
                        EXTENT...,
                    ) 
                c = N_ATTEMPTS
            catch 
                @info "failed $l $ct times"
                ct += 1 
            end
        end

        try 

            tmp = 
                    SimpleSDMPredictor(
                        RasterData(CHELSA2, BioClim),
                        _scenario_to_projection(s);
                        timespan=_scenario_to_year_pair(s),
                        layer=l,
                        EXTENT...,
                    ) 
        catch
            @info "never successed with $l"
        end
    end 
end