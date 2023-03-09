struct Interaction{T<:Site}
    bee::Bee
    plant::Plant
    int::MangalInteraction
    elevation 
    time
end
Base.show(io::IO, int::Interaction{T}) where T = Base.print(io,"🐝 $(int.bee.name) ↔️ 🌷 $(int.plant.name) at $(sitename(T)) on $(monthname(int.time)) $(day(int.time)), $(year(int.time))")

struct BeeData
    bees::Vector{Bee}
    plants::Vector{Plant}
    interactions::Vector{Interaction}
    occurrence::DataFrame
    environment::DataFrame
    cooccurence::DataFrame
end 
Base.show(io::IO, bd::BeeData) = Base.print(io, "Pollination dataset with $(length(bd.interactions)) interactions")

occurrence(bd::BeeData) = bd.occurrence
environment(bd::BeeData) = bd.environment
bees(bd::BeeData) = bd.bees
plants(bd::BeeData) = bd.plants

plant(int::MangalInteraction) = int.to
plant(int::Interaction) = int.plant
plant(bd::BeeData, str) = plants(bd)[findfirst(x->x.name==str, plants(bd))]

bee(int::MangalInteraction) = int.from
bee(int::Interaction) = int.bee
bee(bd::BeeData, str) = bees(bd)[findfirst(x->x.name==str,bees(bd))]
