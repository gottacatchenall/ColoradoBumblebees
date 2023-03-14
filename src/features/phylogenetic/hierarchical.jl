
struct Hierarchical <: Phylogenetic end
outdim(h::Hierarchical, ::Type{Bee}) = 6
outdim(h::Hierarchical, ::Type{Plant}) = 13

function getfeatures(h::Hierarchical, data)
    beestr, plantstr = get_newick()
    beetree, planttree = readnw(beestr), readnw(plantstr)

    species = vcat(bees(data)..., plants(data)...)
    return codes = merge([make_codes(t, species) for t in [beetree, planttree]]...)
end

function make_codes(tree, species)
    height = treeheight(tree)
    dict = Dict()
    for x in Leaves(tree)
        tmp = x
        code = zeros(Float32, height)
        cursor = 1
        while !isnothing(parent(tmp))
            i = findfirst(isequal(tmp), children(parent(tmp))) - 1
            code[end + 1 - cursor] = i
            cursor += 1
            tmp = parent(tmp)
        end

        sp = species[findfirst(s -> x.data.name == s.name, species)]

        merge!(dict, Dict(sp => code))
    end
    return dict
end
