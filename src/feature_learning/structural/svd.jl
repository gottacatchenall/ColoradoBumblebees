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
        x = b[findfirst(x -> x.name == bee_node, b)]
        merge!(dict, Dict(x => L[i, 1:(s.dimensions)]))
    end
    for (i, plant_node) in enumerate(mw.B)
        x = p[findfirst(x -> x.name == plant_node, p)]
        merge!(dict, Dict(x => R[i, 1:(s.dimensions)]))
    end
    return dict
end