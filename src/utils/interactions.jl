interactions(bd::BeeData) = bd.interactions
function interactions(bd::BeeData, sp1, sp2)
    thisbee, thisplant = split(sp1.name, " ")[1] == "Bombus" ? (sp1, sp2) : (sp2, sp1)
    return filter(int -> bee(int) == thisbee && plant(int) == thisplant, interactions(bd))
end
function interactions(bd::BeeData, sp1)
    isbee = split(sp1.name, " ")[1] == "Bombus"

    otherspecies = isbee ? plants(bd) : bees(bd)
    Is = []
    for sp2 in otherspecies
        thisbee = isbee ? sp1 : sp2
        thisplant = isbee ? sp2 : sp1
        Is = vcat(
            Is...,
            findall(
                int -> bee(int) == thisbee && plant(int) == thisplant, interactions(bd)
            )...,
        )
    end
    return interactions(bd)[Is]
end

function metaweb(data::BeeData)
    mangalnet = convert(
        BipartiteNetwork, convert(UnipartiteNetwork, [i.int for i in interactions(data)])
    )

    allbees, allplants = bees(data), plants(data)
    t = [
        allbees[findfirst(x -> x.mangalnode == mangalnet.T[i], allbees)] for
        i in 1:length(mangalnet.T)
    ]
    b = [
        allplants[findfirst(x -> x.mangalnode == mangalnet.B[i], allplants)] for
        i in 1:length(mangalnet.B)
    ]

    return net = BipartiteNetwork(adjacency(mangalnet), t, b)
end
