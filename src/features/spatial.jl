struct KMeansSpatialEmbedding <: Spatial
    k  # embedding dimensions
end

outdim(kmse::KMeansSpatialEmbedding) = 2kmse.k
outdim(kmse::KMeansSpatialEmbedding, ::Union{Type{Bee},Type{Plant}}) = outdim(kmse)

function getfeatures(kmse::KMeansSpatialEmbedding, data)
    occ = occurrence(data)
    species = vcat(bees(data)..., plants(data)...)
    dict = Dict()
    for s in species
        thisspecies = findall(x->x.species==s.name, eachrow(occ))
        merge!(dict, Dict(s=>getfeatures(kmse, occ[thisspecies, :])))
    end
    dict
end

function getfeatures(kmse::KMeansSpatialEmbedding, occ::DataFrame)
    long, lat = occ[!, :longitude], occ[!,:latitude]
    X = hcat(long, lat)'
    if nrow(occ) > kmse.k
        res = kmeans(X, kmse.k)
        return  vcat([res.centers[:,i] for i in 1:kmse.k]...)    
    else
        @info "no work with this species, still dumb hotfix"
        vcat([[mean(occ.longitude), mean(occ.latitude)] for i in 1:kmse.k]...)
    end
end

#=
fig = Figure()
ax1 = Axis(fig[1,1])
scatter!(ax1, long,lat, color=(:blue, 0.15))
scatter!(ax1, [c[1] for c in ctrs], [c[2] for c in ctrs], markersize=20, color=(:orange,0.5))
current_figure()    
=#