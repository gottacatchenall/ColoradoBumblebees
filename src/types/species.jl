abstract type Species end 

struct Bee <: Species
    name::Any
    mangalnode::MangalNode
end
Base.show(io::IO, bee::Bee) = Base.print(io, "🐝 $(bee.name)")

struct Plant <: Species
    name::Any
    mangalnode::MangalNode
end
Base.show(io::IO, plant::Plant) = Base.print(io, "🌷 $(plant.name)")
