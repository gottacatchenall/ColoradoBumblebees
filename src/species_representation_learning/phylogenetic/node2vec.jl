Base.@kwdef struct PhylogeneticNode2Vec <: Phylogenetic
    p = 0.5
    q = 0.5
    number_of_walks = 100
    walk_length = 50
    embedding_dim = 10
end
outdim(pn2v::PhylogeneticNode2Vec, ::Union{Type{Bee},Type{Plant}}) = pn2v.embedding_dim
outdim(pn2v::PhylogeneticNode2Vec) = pn2v.embedding_dim

function _embed(data::BeeData, pn2v::PhylogeneticNode2Vec)
    beetree, planttree = readnw.(load_newick())
    return merge([_run_node2vec(pn2v, t, data) for t in [beetree, planttree]]...)
end

function _run_node2vec(pn2v::PhylogeneticNode2Vec, tree, data)
    g, Ileaf, species_names = tree_to_graph(tree)
    walks = simulate_walks(
        g; num_walks=pn2v.number_of_walks, p=pn2v.p, q=pn2v.q, walk_length=pn2v.walk_length
    )
    vecs = learn_embeddings(walks; outdim=outdim(pn2v))
    vecs = vecs[:, Ileaf]

    dict = Dict()
    for (i, sp) in enumerate(species_names)
        isbee = split(sp, " ")[1] == "Bombus"
        f = isbee ? bee : plant
        speciesobj = f(data, sp)
        merge!(dict, Dict(speciesobj => vecs[:, i]))
    end
    return dict
end

function tree_to_graph(tree)
    g = SimpleGraph()
    dict = Dict()

    leaf_indices = []
    leaf_species_names = []

    for (i, n) in enumerate(PreOrderDFS(tree))
        NewickTree.isleaf(n) && push!(leaf_indices, i)
        a = split(n.data.name, "_")
        NewickTree.isleaf(n) && push!(leaf_species_names, string(a[1], " ", a[2]))
        add_vertex!(g)
        merge!(dict, Dict(n.id => i))
    end
    for (i, n) in enumerate(PreOrderDFS(tree))
        for c in children(n)
            cid = c.id
            j = dict[cid]
            add_edge!(g, i, j)
        end
    end
    return g, leaf_indices, leaf_species_names
end
