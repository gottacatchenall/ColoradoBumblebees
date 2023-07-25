# Addapted from github.com/ollin18/Node2Vec.jl to work with Graphs.jl
function simulate_walks(g; num_walks=50, kwargs...)
    walks = Array{Array}(undef, 0)
    nodes = collect(vertices(g))
    for _ in 1:num_walks
        nodes = Random.shuffle(nodes)
        for node in nodes
            push!(walks, node2vec(g, node; kwargs...))
        end
    end
    return walks
end

function learn_embeddings(walks; outdim::Int=100)
    str_walks = map(x -> string.(x), walks)

    # THIS NEEDS TO BE THREAD SAFE
    rpath = "/tmp"
    the_walks = joinpath(rpath, "str_walk_THREAD$(Threads.threadid()).txt")
    the_vecs = joinpath(rpath, "str_walk-vec_THREAD$(Threads.threadid()).txt")

    writedlm(the_walks, str_walks)
    word2vec(the_walks, the_vecs; verbose=true, size=outdim)
    model = wordvectors(the_vecs)
    return model.vectors
end

function node2vec(g, start; p=0.5, q=0.5, walk_length=10)
    no, ed = preprocess_transition_probs(g, p, q)
    walk = [start]
    while length(walk) < walk_length
        current = last(walk)
        neigh = neighbors(g, current)
        if length(neigh) > 0
            if length(walk) == 1
                append!(walk, neigh[alias_draw(no[current][1], no[current][2])])
            else
                prev = walk[length(walk) - 1]
                next = neigh[alias_draw(ed[(prev, current)][1], ed[(prev, current)][2])]
                append!(walk, next)
            end
        else
            break
        end
    end
    return walk
end

function preprocess_transition_probs(g, p, q)
    alias_nodes = Dict()
    for node in vertices(g)
        neigh = neighbors(g, node)
        alias_nodes[node] = alias_setup([1.0 / length(neigh) for n in neigh])
    end
    alias_edges = Dict()
    for e in edges(g)
        alias_edges[(e.src, e.dst)] = get_alias_edge(g, e.src, e.dst, p, q)
        alias_edges[(e.dst, e.src)] = get_alias_edge(g, e.dst, e.src, p, q)
    end
    return alias_nodes, alias_edges
end

function get_alias_edge(g, src, dst, p, q)
    unnormalized_probs = Array{Float64}(undef, 0)
    for dst_nbr in neighbors(g, dst)
        if dst_nbr == src
            append!(unnormalized_probs, 1 / p)
        elseif has_edge(g, dst_nbr, src)
            append!(unnormalized_probs, 1)
        else
            append!(unnormalized_probs, 1 / q)
        end
    end
    norm_const = sum(unnormalized_probs)
    normalized_probs = unnormalized_probs / norm_const
    return alias_setup(normalized_probs)
end

function alias_setup(probs)
    K = length(probs)
    q = zeros(K)
    J = zeros(Int, K)

    smaller = Array{Float64}(undef, 0)
    larger = Array{Float64}(undef, 0)

    for (kk, prob) in enumerate(probs)
        q[kk] = K * prob
        if q[kk] < 1
            push!(smaller, kk)
        else
            push!(larger, kk)
        end
    end

    while length(smaller) > 0 && length(larger) > 0
        small = Int(pop!(smaller))
        large = Int(pop!(larger))

        J[small] = large
        q[large] = q[large] + q[small] - 1.0
        if q[large] < 1.0
            append!(smaller, large)
        else
            append!(larger, large)
        end
    end
    return J, q
end

function alias_draw(J, q)
    K = length(J)
    kk = Int(ceil(rand() * K))
    if rand() < q[kk]
        return kk
    else
        return J[kk]
    end
end
