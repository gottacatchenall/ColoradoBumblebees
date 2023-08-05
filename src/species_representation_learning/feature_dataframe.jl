
function feature_dataframe(data, rep::SpeciesRepresentations)
    label_df = label_dataframe(data)
    
    feature_dims = 2outdim(rep.embed_model)
    emb = rep.embedding

    feature_df = _initialize_df(rep.embed_model, feature_dims, label_df)

    for (ri, r) in enumerate(eachrow(label_df))
        b, p = r.bee, r.plant
        beevec, plantvec = emb[b], emb[p]
        for i in eachindex(beevec)
            feature_df[ri, "$(string(typeof(rep.embed_model)))_BEE$i"] = beevec[i]
        end
        for i in eachindex(plantvec)
            feature_df[ri, "$(string(typeof(rep.embed_model)))_PLANT$i"] = plantvec[i]
        end
    end

    feature_df = MLJ.transform(fit!(machine(Standardizer(), feature_df)),feature_df)
    return hcat(label_df, feature_df)
end

function feature_dataframe(data, reps::Vector{S}) where S<:SpeciesRepresentations
    label_df = label_dataframe(data)
    feature_df = hcat([_initialize_df(r.embed_model, 2outdim(r.embed_model), label_df) for r in reps]...)
    for (ri, r) in enumerate(eachrow(label_df))
        b, p = r.bee, r.plant
        for rep in reps
            emb = rep.embedding
            beevec, plantvec = emb[b], emb[p]
            for i in eachindex(beevec)
                feature_df[ri, "$(string(typeof(rep.embed_model)))_BEE$i"] = beevec[i]
            end
            for i in eachindex(plantvec)
                feature_df[ri, "$(string(typeof(rep.embed_model)))_PLANT$i"] = plantvec[i]
            end
        end 
    end
    feature_df = MLJ.transform(fit!(machine(Standardizer(), feature_df)),feature_df)
    hcat(label_df, feature_df)
end

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

function _initialize_df(emb_model, feature_dims, label_df)
    per_species_dim = feature_dims ÷ 2
    
    cols = []
    push!(cols, ["$(string(typeof(emb_model)))_BEE$fi" => zeros(nrow(label_df)) for fi in 1:per_species_dim]...)
    push!(cols, ["$(string(typeof(emb_model)))_PLANT$fi" => zeros(nrow(label_df)) for fi in 1:per_species_dim]...)
    DataFrame(cols...)
end 

