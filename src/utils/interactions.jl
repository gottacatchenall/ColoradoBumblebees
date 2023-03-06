interactions(bd::BeeData) = bd.interactions
function interactions(bd::BeeData, sp1, sp2)
    thisbee, thisplant = split(sp1.name, " ")[1] == "Bombus" ? (sp1, sp2) : (sp2, sp1)
    filter(int ->  bee(int) == thisbee && plant(int) == thisplant, interactions(bd))
end
function interactions(bd::BeeData, sp1)
    isbee = split(sp1.name, " ")[1] == "Bombus"

    otherspecies = isbee ? plants(bd) : bees(bd)
    Is = []
    for sp2 in otherspecies
        thisbee = isbee ? sp1 : sp2
        thisplant = isbee ? sp2 : sp1
        Is = vcat(Is..., findall(int ->  bee(int) == thisbee && plant(int) == thisplant, interactions(bd))...)
    end
    interactions(bd)[Is]
end