using DrWatson
@quickactivate :ColoradoBumblebees

data = load_data()

edgelist = label_dataframe(data)

allbees, allplants = unique(bees(data)), unique(plants(data))

allspecies = vcat(allbees, allplants)

id_df = DataFrame(species=[s.name for s in allspecies], node_id=[i for i in eachindex(allspecies)])

# TODO need each id to be unique. Two dicts: edgelist and translating species_id
# to speices name 
edgelist_df = DataFrame()
edgelist_df.bee = [findfirst(isequal(s), allspecies) for s in edgelist.bee]
edgelist_df.plant = [findfirst(isequal(s), allspecies) for s in edgelist.plant]


vgae_dir = joinpath(artifactdir(), "vgae")
mkpath(vgae_dir)

CSV.write(joinpath(vgae_dir, "edgelist.csv"), edgelist_df)
CSV.write(joinpath(vgae_dir, "node_ids.csv"), id_df)

