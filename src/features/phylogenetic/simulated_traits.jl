
function simtraits(tr; numsims = 1000)
    df = DataFrame(nodename = getnodename.(tr, traversal(tr, preorder)))
    for rep = 1:numsims
        σ = rand(Exponential(1.0))
        rand!(BrownianTrait(tr, "Trait($rep)"; σ = σ), tr)
        d = DataFrame(
            nodename = getnodename.(tr, traversal(tr, preorder)),
            trait = getnodedata.(tr, traversal(tr, preorder), "Trait($rep)"),
        )
        df[!, "trait$(rep)"] = d.trait
    end
    df

    leafnames = collect(nodenamefilter(isleaf, tr))
    filter(row -> row.nodename in leafnames, df)
end

function makepca(tree; kwargs...)
    df = simtraits(tree; kwargs...)
    ntaxa, ntraits = nrow(df), ncol(df) - 1
    pcainput = zeros(ntaxa, ntraits)

    for (i, row) in enumerate(eachrow(df))
        for colindex = 2:ncol(df)
            pcainput[i, colindex-1] = row[colindex]
        end
    end

    df, pcainput
end

function beephylogenyfeatures(; phylodims = 5, kwargs...)
    beetree = open(parsenewick, joinpath(PROJECT_ROOT, "data", "artifacts", "bee_tree.tre"))
    df, pcainput = makepca(beetree; kwargs...)

    df = _replace_bee_codenames(df)
    pca = fit(PCA, pcainput')
    m = MultivariateStats.transform(pca, pcainput')
    pcadf = DataFrame(species = df.species)

    outdim = size(m, 1)
    lastdim = min(phylodims, outdim)
    for i = 1:lastdim
        pcadf[!, Symbol("dim$i")] = m[i, :]
    end
    return pcadf
end
