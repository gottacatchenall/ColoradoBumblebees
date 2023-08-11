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



# This is essential for saved and loaded data artifacts to read as == to objects
# retrieved via `bee(data, species_name)`

Base.isequal(a::A, b::B) where {A<:Species,B<:Species} = a.name == b.name
Base.:(==)(a::A, b::B) where {A<:Species,B<:Species} = a.name == b.name 

Base.haskey(d::Dict{S,V}, key::S) where {S<:Species, V} = any(isequal(key), keys(d))

# This makes it linear. Doesn't really matter, but kind of annyoing
Base.getindex(d::Dict{S,V}, key::S) where {S<:Species, V} = begin
    k,v = collect(keys(d)), collect(values(d))
    v[findfirst(isequal(key), k)]
end