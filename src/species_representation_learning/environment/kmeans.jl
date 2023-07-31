Base.@kwdef struct KMeansEnvironmentEmbedding <: Environment
    k::Any = 4 # embedding dimensions
end

outdim(kmee::KMeansEnvironmentEmbedding) = 3kmee.k
outdim(kmee::KMeansEnvironmentEmbedding, ::Union{Type{Bee},Type{Plant}}) = outdim(kmee)

function ColoradoBumblebees.getfeatures(kmee::KMeansEnvironmentEmbedding, data)
    env = environment(data)
    species = vcat(bees(data)..., plants(data)...)
    dict = Dict()
    for s in species
        thisspecies = findall(x -> x.species == s.name, eachrow(env))
        merge!(dict, Dict(s => getfeatures(kmee, env[thisspecies, :])))
    end

    ext = vcat(values(dict)...)
    m, M = extrema(ext)
    for (k, v) in dict
        dict[k] = (v .- m) ./ (M - m)
    end

    return dict
end

function getfeatures(kmse::KMeansEnvironmentEmbedding, env::DataFrame)
    X = Matrix(env[!, 2:end])'
    #X = hcat(long, lat)'
    if nrow(env) > kmse.k
        res = kmeans(X, kmse.k)
        return vcat([res.centers[:, i] for i in 1:(kmse.k)]...)
    else
        @info "no work with this species, still dumb hotfix"
        return vcat([[mean(env[!, "w_BIO_$i"]) for i in 1:size(X, 1)] for i in 1:(kmse.k)]...)
        #vcat([[mean(occ.longitude), mean(occ.latitude)] for i in 1:kmse.k]...)
    end
end
