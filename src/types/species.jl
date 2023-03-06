struct Bee 
    name
    mangalnode::MangalNode
end
Base.show(io::IO, bee::Bee) = Base.show(io, "🐝 $(bee.name)")

struct Plant
    name
    mangalnode::MangalNode
end 
Base.show(io::IO, plant::Plant) = Base.show(io, "🌷 $(plant.name)")
