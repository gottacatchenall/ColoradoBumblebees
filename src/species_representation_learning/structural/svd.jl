Base.@kwdef struct MetawebSVD <: Structural
    truncation_dims = 8
    embed_dims = 8
end
outdim(svd::MetawebSVD) = svd.embed_dims
outdim(svd::MetawebSVD, ::Union{Type{Bee},Type{Plant}}) = outdim(svd)

function _embed(data::BeeData, s::MetawebSVD)
    trunc_dims =  s.truncation_dims
    emb_dims = s.embed_dims
    mw = metaweb(data)
    truncated_M = _get_truncated_metaweb(mw, trunc_dims)

    bee_emb, _, plant_emb = svd(truncated_M)
    bee_names, plant_names = mw.T, mw.B

    embedding_dict = merge(Dict(zip(bee_names, eachrow(bee_emb))), Dict(zip(plant_names, eachrow(plant_emb))))
    _create_embedding_dict(data, embedding_dict, emb_dims)
end 
