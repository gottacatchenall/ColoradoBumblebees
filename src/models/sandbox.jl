function feature_dataframe(data, embeddings)
    labeldf = label_dataframe(data)

    per_species_dim = sum(outdim.(embeddings))

    featdf = DataFrame([["BEE$i"=>[] for i in 1:per_species_dim]...,["PLANT$i"=>[] for i in 1:per_species_dim]... ])
    
    feats = [getfeatures(e, data) for e in embeddings]

    for r in eachrow(labeldf)
        b,p = r.bee, r.plant
        cursor = 1
        for f in feats
            beevec, plantvec = f[b], f[p]
            for i in 1:length(f[b])
                push!(featdf[!, "BEE$cursor"], beevec[i])
                push!(featdf[!, "PLANT$cursor"], plantvec[i])
                cursor += 1
            end
        end
    end
    hcat(labeldf, featdf)
end


# each row is pair of species
function label_dataframe(data)
    allpairs = [(b,p) for b in bees(data), p in plants(data)]

    df = DataFrame(bee=[], plant=[], interaction=[])
    for I in eachindex(allpairs)
        b,p = allpairs[I]
        push!(df.bee, b)
        push!(df.plant,p)
        push!(df.interaction, length(interactions(data,b,p)) > 0 )
    end
    df
end