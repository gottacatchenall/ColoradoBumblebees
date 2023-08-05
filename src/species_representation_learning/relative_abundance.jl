struct Pooled <: RelativeAbundance end

outdim(p::Pooled) = 1

function _embed(data::BeeData, ::Pooled)
    ints = interactions(data)

    b, p = bees(data), plants(data)

    allints = length(ints)

    beeabd = Dict([
        (thisbee => [length(interactions(data, thisbee)) / allints]) for thisbee in b
    ])
    plantabd = Dict([
        (thisplant => [length(interactions(data, thisplant)) / allints]) for thisplant in p
    ])

    return merge(beeabd, plantabd)
end
