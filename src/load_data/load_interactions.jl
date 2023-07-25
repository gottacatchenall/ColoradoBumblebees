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
    b,p = bees(data), plants(data)
    ints = ColoradoBumblebees.interactions(data)
    adj = zeros(Bool, length(b), length(p))
    for int in ints
        Ib = findfirst(isequal(int.bee), b) 
        Ip = findfirst(isequal(int.plant), p) 
        adj[Ib,Ip] = 1
    end
    BipartiteNetwork(adj, [bi.name for bi in b], [i.name for i in p])
end
