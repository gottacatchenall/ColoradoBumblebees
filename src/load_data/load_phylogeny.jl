function load_newick()
    beestr = read(datadir("public", "phylogeny", "newick_trees", "bee_tree.newick"), String)
    plantstr = read(datadir("public", "phylogeny", "newick_trees", "plant_tree.newick"), String)

    bee_df = CSV.read(datadir("public", "phylogeny", "raw_sequences", "bee_sequences.csv"), DataFrame)
    rename!(bee_df, :Bee => :species)
    plant_df = CSV.read(datadir("public", "phylogeny", "raw_sequences", "plant_sequences.csv"), DataFrame)
    rename!(plant_df, :Plants => :species)
    
    codenames = Dict([
        r["Species Codename"] => r["species"] for r in eachrow(
            vcat(
                bee_df[!, ["species", "Species Codename"]],
                plant_df[!, ["species", "Species Codename"]],
            ),
        )
    ])

    for (k, v) in codenames
        spl = split(v, " ")
        str = string(spl[1], "_", spl[2])

        plantstr = replace(plantstr, "$k:" => "$str:")
        beestr = replace(beestr, "$k:" => "$str:")
    end
    return beestr, plantstr
end
