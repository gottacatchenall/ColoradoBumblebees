using DrWatson
@quickactivate :ColoradoBumblebees

data = load_data()

dirpaths = joinpath(artifactdir(), "vgae", "fits")
filepaths = readdir(dirpaths)
 
a = joinpath.(artifactdir(), "vgae", "fits", filepaths)

embed_dfs = CSV.read.(a, DataFrame)


input_features = [parse(Int32, split(x,"_")[3]) for x in filepaths]
embed_dims = [parse(Int32, split(split(x,"_")[end], ".")[begin]) for x in filepaths]

models = [contains(x, "VGAE") ? GraphAutoencoder{Variational}(input_features[i], embed_dims[i]) : GraphAutoencoder{Standard}(input_features[i], embed_dims[i]) for (i,x) in enumerate(filepaths)]

reps = []

for (i,df) in enumerate(embed_dfs)
    dict = Dict{Species, Vector}()
    for r in eachrow(df)
        sp = contains(r.species, "Bombus") ? bee(data, r.species) : plant(data, r.species)
        merge!(dict, Dict(sp => [x[1] for x in r[2:end-1]]))
    end
    push!(reps, SpeciesRepresentations(models[i], dict, nothing))
end 


for r in reps
    ColoradoBumblebees.save(r)
end