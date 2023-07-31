Base.@kwdef struct LFSVD <: Structural
    α = [0.25 0.25 0.25 0.25]
    truncation_dims = 8
    embed_dims = 8
end
outdim(svd::LFSVD) = svd.embed_dims
outdim(svd::LFSVD, ::Union{Type{Bee},Type{Plant}}) = outdim(svd)

function getfeatures(s::LFSVD, data)
    trunc_dims =  s.truncation_dims
    emb_dims = s.embed_dims

    mw = ColoradoBumblebees.metaweb(data)
    lf = linearfilter(mw, s.α)

    truncated_M = _get_truncated_metaweb(lf, trunc_dims)


    bee_emb, _, plant_emb = svd(truncated_M)
    embedding_dict = merge(Dict(zip(mw.T, eachrow(bee_emb))), Dict(zip(mw.B, eachrow(plant_emb))))
    _create_embedding_dict(data, embedding_dict, emb_dims)
end

function linearfilter(N, α) 
    α = α ./ sum(α)

    gen = [mean(x) for x in eachrow(N.edges)]
    vul =  [mean(x) for x in eachcol(N.edges)]
    co = mean(N.edges)

    P = zeros(size(N))
    for i in axes(N.edges, 1)
        for j in axes(N.edges, 2)
            P[i, j] = α[1] * N[i, j] + α[2] * gen[i] + α[3] * vul[j] + α[4] * co
        end
    end

    clamp!(P, zero(typeof(P[begin])), one(typeof(P[begin])))

    return BipartiteProbabilisticNetwork(P, N.T, N.B)
end