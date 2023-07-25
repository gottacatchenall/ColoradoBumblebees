using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

using CSV, DataFrames

data = load_data()

phen = load_phenology(data)

df = DataFrame(species=[], time=[], abundance=[])

for (k,v) in phen
    species_name = k.name
    for (t,x) in enumerate(v) 
        push!(df.species, species_name)
        push!(df.time, t)
        push!(df.abundance, x)
    end
end

df

CSV.write(scriptsdir("supplement", "temporal_autoencoder_model_tuning", "pheno.csv"), df)
