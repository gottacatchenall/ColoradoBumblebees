struct SpeciesRepresentations{T,P<:Union{String, Nothing}}
    embed_model::T
    embedding::Dict{Species, Vector}
    path::P
end

Base.show(io::IO, sr::SpeciesRepresentations{T,P}) where {T,P} = Base.print(io, "📊 $(supertype(T)) species representations using $(typeof(sr.embed_model))")


function representations(data, model)
    feat_dict = embed(data, model)
    SpeciesRepresentations(model, feat_dict, nothing)
end

path(sr::SpeciesRepresentations) = joinpath(
    artifactdir(), 
    "species_representations", 
    "SpeciesRepresenetation_$(typeof(sr.embed_model))_"*savename(sr.embed_model)
)


DrWatson.default_prefix(sr::SpeciesRepresentations{T,P}) where {T,P} = "SpeciesRepresentation_"
DrWatson.savename(sr::SpeciesRepresentations{T,P}) where {T,P} = string(T)*"("*savename(sr.embed_model)*")"

