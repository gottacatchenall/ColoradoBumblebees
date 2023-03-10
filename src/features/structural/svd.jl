Base.@kwdef struct MetawebSVD <: Structural
    dimensions = 8
end
outdim(svd::MetawebSVD) = svd.dimensions
outdim(svd::MetawebSVD, ::Union{Type{Bee},Type{Plant}}) = outdim(svd)

function getfeatures(s::MetawebSVD, data)
    mw = ColoradoBumblebees.metaweb(data)
    L, Λ, R = svd(adjacency(mw))

    b, p = bees(data), plants(data)

    dict = Dict()
    for (i, bee_node) in enumerate(mw.T)
        x = b[findfirst(x -> x.name == bee_node.name, b)]
        merge!(dict, Dict(x => L[i, 1:(s.dimensions)]))
    end
    for (i, plant_node) in enumerate(mw.B)
        x = p[findfirst(x -> x.name == plant_node.name, p)]
        merge!(dict, Dict(x => R[i, 1:(s.dimensions)]))
    end
    return dict
end

Base.@kwdef struct CooccurencePCA <: Structural
    dimensions = 16
end
outdim(pca::CooccurencePCA) = pca.dimensions

function getfeatures(pca::CooccurencePCA, data)
    C, proj = pca_coocc(data)

    sp = vcat(bees(data)..., plants(data)...)

    dict = Dict()

    for s in sp
        col = findall(x -> x == s.name, C.S)[begin]
        merge!(dict, Dict(s => proj[1:(pca.dimensions), col]))
    end
    return dict
end

function pca_coocc(data)
    C = cooccurence(data)

    B = BipartiteNetwork(C, [b.name for b in bees(data)], [p.name for p in plants(data)])
    N = convert(UnipartiteNetwork, B)
    K = EcologicalNetworks.mirror(N)
    pc = MultivariateStats.fit(PPCA, Float64.(Array(K.edges)))
    pr = MultivariateStats.transform(pc, Float64.(Array(K.edges)))
    return K, pr
end

#=
K, proj = pca_coocc()
using CairoMakie
x,y = [proj[2,i] for i in 1:size(proj,2)], [proj[3,i] for i in 1:size(proj,2)]
scatter(x,y) =#
