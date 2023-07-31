Base.@kwdef struct SimulatedTraits <: Phylogenetic
    numtraits = 1000
    variance_distribution = Exponential(1.0)
    truncated_dims = 8
end
outdim(st::SimulatedTraits, ::Union{Type{Bee},Type{Plant}}) = st.truncated_dims

function getfeatures(st::SimulatedTraits, data)
    dict = Dict()

    species = vcat(bees(data), plants(data))

    bee_tree, plant_tree = parsenewick.(load_newick())
    df = vcat(
        _simulated_trait_features(st, bee_tree), _simulated_trait_features(st, plant_tree)
    )
    
    for r in eachrow(df)
        Ifirst = findfirst(s -> s.name == r.species, species)
        !isnothing(Ifirst) && merge!(dict, Dict(species[Ifirst] => [x for x in r[2:end]]))
    end
    return dict
end

function _simulated_trait_features(st, tree; kwargs...)
    truncated_dims = st.truncated_dims

    df, pcainput = _make_pca_input(tree; kwargs...)
    spl = split.(df.nodename, "_")

    species_fixed = [string(s[1], " ", s[2]) for s in spl]

    pca = MultivariateStats.fit(PCA, pcainput')
    m = MultivariateStats.transform(pca, pcainput')
    pcadf = DataFrame(; species=species_fixed)
    outdim = size(m, 1)
    lastdim = min(truncated_dims, outdim)
    for i in 1:lastdim
        pcadf[!, Symbol("dim$i")] = m[i, :]
    end
    return pcadf
end

function _simulate_traits(tree; numsims=1000)
    df = DataFrame(; nodename=getnodename.(tree, traversal(tree, preorder)))
    for rep in 1:numsims
        σ² = rand(Exponential(1.0))
        rand!(BrownianTrait(tree, "Trait($rep)"; σ²=σ²), tree)
        d = DataFrame(;
            nodename=getnodename.(tree, traversal(tree, preorder)),
            trait=getnodedata.(tree, traversal(tree, preorder), "Trait($rep)"),
        )
        df[!, "trait$(rep)"] = d.trait
    end
    df

    leafnames = collect(nodenamefilter(Phylo.isleaf, tree))
    return filter(row -> row.nodename in leafnames, df)
end

function _make_pca_input(tree; kwargs...)
    df = _simulate_traits(tree; kwargs...)
    ntaxa, ntraits = nrow(df), ncol(df) - 1
    pcainput = zeros(ntaxa, ntraits)

    for (i, row) in enumerate(eachrow(df))
        for colindex in 2:ncol(df)
            pcainput[i, colindex - 1] = row[colindex]
        end
    end
    return df, pcainput
end
