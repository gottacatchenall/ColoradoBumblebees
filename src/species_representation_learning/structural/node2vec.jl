@Base.kwdef struct MetawebNode2Vec <: Structural
    p = 0.5
    q = 0.5
    number_of_walks = 100
    walk_length = 50
    embedding_dim = 10
end

outdim(mn2v::MetawebNode2Vec, ::Union{Type{Bee},Type{Plant}}) = mn2v.embedding_dim
outdim(mn2v::MetawebNode2Vec) = mn2v.embedding_dim

function _embed(data::BeeData, mn2v::MetawebNode2Vec)

    M = EcologicalNetworks.mirror(convert(UnipartiteNetwork, metaweb(data)))


    adj, spnames = adjacency(M), M.S

    g = SimpleGraph(adj)

    walks = ColoradoBumblebees.simulate_walks(
        g; 
        num_walks=mn2v.number_of_walks, 
        p=mn2v.p, 
        q=mn2v.q, 
        walk_length=mn2v.walk_length
    )

    vecs = ColoradoBumblebees.learn_embeddings(walks; outdim=outdim(mn2v))

    dict = Dict()
    for (i, sp) in enumerate(spnames)
        isbee = split(sp, " ")[1] == "Bombus"
        f = isbee ? bee : plant
        speciesobj = f(data, sp)
        !isnothing(speciesobj) && merge!(dict, Dict(speciesobj => vecs[:, i]))
    end
    return dict
end
