const SR = Union{SpeciesRepresentations,Vector{S}} where S<:SpeciesRepresentations

struct ClassificationFit{M<:ClassificationModel,V<:SR}
    model::M
    representation::V
    predictions::DataFrame  # todo: add column that is whether int was in train or test
    fit_stats::Dict
end
Base.show(io::IO, cf::ClassificationFit{M,V}) where {M,V} = Base.print(io, "ClassificationFit with $M on $(_name_reptypes(cf.representation))")

struct BatchFit{M,V}
    fits::Vector{ClassificationFit{M,V}}
end

_name_reptypes(::ClassificationFit{M,S}) where {M,S} = typeof(S)
_rep_type(::SpeciesRepresentations{T}) where T = T
_name_reptypes(cf::ClassificationFit{M,Vector{S}}) where {M,S} = string("(", ["$(_rep_type(i)), " for i in cf.representation]..., ")")
Base.show(io::IO, bf::BatchFit{M,V}) where {M,V} = Base.print(io, "BatchFit with $(length(bf.fits)) replicates fit with $M on $(_name_reptypes(bf.fits[1]))")


ColoradoBumblebees.path(bf::BatchFit) = joinpath(
    artifactdir(), 
    "classification_fits", 
    savename(bf)
)


function DrWatson.savename(bf::BatchFit{M,Vector{S}}) where {M,S}
    str = savename(bf.fits[1].model)

    str *= "}_Embedding{"
    for r in bf.fits[1].representation[1:end-1]
        str *= savename(r)*"_"
    end
    str *= savename(bf.fits[1].representation[end]) * "}"
    str
end 


function DrWatson.savename(bf::BatchFit{M,S}) where {M,S}
    str = savename(bf.fits[1].model)
    str *= "}_Embedding{"*savename(bf.fits[1].representation)*"}"
    str
end 

 