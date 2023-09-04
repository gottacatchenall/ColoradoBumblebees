const SR = Union{SpeciesRepresentations,Vector{S}} where S<:SpeciesRepresentations

struct ClassificationFit{M<:ClassificationModel,V<:SR}
    model::M
    representation::V
    predictions::DataFrame  # todo: add column that is whether int was in train or test
    fit_stats::Dict
end
Base.show(io::IO, cf::ClassificationFit{M,V}) where {M,V} = Base.print(io, "ClassificationFit with $M on $(representation(cf))")

model(cf::ClassificationFit) = cf.model
fit_stats(cf::ClassificationFit) = cf.fit_stats
predictions(cf::ClassificationFit) = cf.predictions
representation(cf::ClassificationFit) = cf.representation
prauc(cf::ClassificationFit) = haskey(fit_stats(cf), "prauc") ? fit_stats(cf)["prauc"] : fit_stats(cf)[:prauc]
rocauc(cf::ClassificationFit) = haskey(fit_stats(cf), "rocauc") ? fit_stats(cf)["rocauc"] : fit_stats(cf)[:rocauc]

struct BatchFit{M,V}
    fits::Vector{ClassificationFit{M,V}}
end

model(bf::BatchFit) = model(bf.fits[1])
fit_stats(bf::BatchFit) = fit_stats.(bf.fits)
predictions(bf::BatchFit) = predictions.(bf.fits)
representation(bf::BatchFit) = representation(bf.fits[1])
praucs(bf::BatchFit) = prauc.(bf.fits)
rocaucs(bf::BatchFit) = rocauc.(bf.fits)


_name_reptypes(::ClassificationFit{M,S}) where {M,S} = typeof(S)
_rep_type(::SpeciesRepresentations{T}) where T = T
_name_reptypes(cf::ClassificationFit{M,Vector{S}}) where {M,S} = string("(", ["$(_rep_type(i)), " for i in cf.representation]..., ")")
Base.show(io::IO, bf::BatchFit{M,V}) where {M,V} = Base.print(io, "BatchFit with $(length(bf.fits)) replicates fit with $M on $(_name_reptypes(bf.fits[1]))")


function ColoradoBumblebees.path(bf::BatchFit{M,V}) where {M,V}
    subdir = V <: Vector ? "multiple_representations" :  "single_representation" 

    joinpath(
        artifactdir(), 
        "classification_fits", 
        subdir, 
        savename(bf)
    )
end


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

 