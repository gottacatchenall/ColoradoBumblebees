Base.@kwdef struct LFSVD <: Structural
    α = [0.25 0.25 0.25 0.25]
    dimensions = 8
end
outdim(svd::LFSVD) = svd.dimensions
outdim(svd::LFSVD, ::Union{Type{Bee},Type{Plant}}) = outdim(svd)

function getfeatures(s::LFSVD, data)
    mw = ColoradoBumblebees.metaweb(data)

    lf = linearfilter(mw, s.α)

    L, Λ, R = svd(adjacency(lf))

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