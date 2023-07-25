function feature_dataframe(data, embeddings)
    labeldf = label_dataframe(data)

    feats = [getfeatures(e, data) for e in embeddings]

    cols = []

    for i in eachindex(embeddings)
        e = embeddings[i]
        per_bee_dim = outdim(e, Bee)
        per_plant_dim = outdim(e, Plant)
        for fi in per_bee_dim
            @info "$(string(typeof(e)))_BEE$fi"
        end
        push!(
            cols,
            [
                "$(string(typeof(e)))_BEE$fi" => zeros(nrow(labeldf)) for fi in 1:per_bee_dim
            ]...,
        )
        push!(
            cols,
            [
                "$(string(typeof(e)))_PLANT$fi" => zeros(nrow(labeldf)) for
                fi in 1:per_plant_dim
            ]...,
        )
    end

    featdf = DataFrame(cols...)

    for (ri, r) in enumerate(eachrow(labeldf))
        b, p = r.bee, r.plant
        for (fi, f) in enumerate(feats)
            beevec, plantvec = f[b], f[p]
            for i in 1:length(f[b])
                featdf[ri, "$(string(typeof(embeddings[fi])))_BEE$i"] = beevec[i]
            end
            for i in 1:length(f[p])
                featdf[ri, "$(string(typeof(embeddings[fi])))_PLANT$i"] = plantvec[i]
            end
        end
    end
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
