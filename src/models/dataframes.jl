function feature_dataframe(data, embeddings)
    labeldf = label_dataframe(data)

    per_bee_dim = sum([outdim(e, Bee) for e in embeddings])
    per_plant_dim = sum([outdim(e, Plant) for e in embeddings])

    featdf = DataFrame([
        ["BEE$i" => zeros(Float32, nrow(labeldf)) for i in 1:per_bee_dim]...,
        ["PLANT$i" => zeros(Float32, nrow(labeldf)) for i in 1:per_plant_dim]...,
    ])
    feats = [getfeatures(e, data) for e in embeddings]

    @info per_bee_dim, per_plant_dim

    for (ri, r) in enumerate(eachrow(labeldf))
        b, p = r.bee, r.plant
        for f in feats
            beevec, plantvec = f[b], f[p]
            for i in 1:length(f[b])
                featdf[ri, "BEE$i"] = beevec[i]
            end
            for i in 1:length(f[p])
                featdf[ri, "PLANT$i"] = plantvec[i]
            end
        end
    end

    @info nrow(featdf), nrow(labeldf)

    return hcat(labeldf, featdf)
end

# each row is pair of species
function label_dataframe(data)
    allpairs = [(b, p) for b in bees(data), p in plants(data)]

    df = DataFrame(; bee=[], plant=[], interaction=[])
    for I in eachindex(allpairs)
        b, p = allpairs[I]
        push!(df.bee, b)
        push!(df.plant, p)
        push!(df.interaction, length(interactions(data, b, p)) > 0)
    end
    return df
end
