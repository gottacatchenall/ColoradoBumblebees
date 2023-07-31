# embedding_dict is Dict{String, Vector}
function _create_embedding_dict(data, embedding_dict, emb_dims)
    dict = Dict()

    bee_nodes, plant_nodes = bees(data), plants(data)
    for (species_name, embed) in embedding_dict
        sp = contains(species_name, "Bombus") ? 
            bee_nodes[findfirst(x -> x.name == species_name, bee_nodes)] : 
            plant_nodes[findfirst(x -> x.name == species_name, plant_nodes)] 
        merge!(dict, Dict(sp => embed[1:emb_dims]))
    end

    return dict
end 


function _get_truncated_metaweb(mw, trunc_dims)
    M = adjacency(mw)
    L, Λ, R = svd(M)
    truncated_Λ = diagm(vcat(Λ[1:trunc_dims], [0 for _ in 1:(length(Λ)-trunc_dims)]...))
    truncated_M = L*truncated_Λ*R'
end 

