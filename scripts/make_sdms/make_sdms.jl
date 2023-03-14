using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

make_sdms(data, GaussianBRT(); cluster=true)


#=
pl = plant(data, "Hypochaeris radicata")

scenario = scenarios()[2]

cl = SimpleSDMPredictor(
            RasterData(CHELSA2, BioClim),
            _scenario_to_projection(scenario);
            timespan=_scenario_to_year_pair(scenario),
            layer="BIO1",
            EXTENT...,
        )


layer = similar(convert(Float32, cl))
layer.grid .= 0
convert_occurrence_to_tif!(pl, layer)
pres, abs = get_pres_and_abs(convert(Bool, layer))


get_features_and_labels(pres, abs, [layer])


findall(!isnothing, abs.grid)

allspecies = vcat(bees(data)...,plants(data)...)
for s in allspecies
    convert_occurrence_to_tif!(s, layer)
    @info s, sum(layer.grid)
end 

=#