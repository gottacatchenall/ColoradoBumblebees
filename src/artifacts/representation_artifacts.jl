function ColoradoBumblebees.save(x::SpeciesRepresentations{T}) where T
    outdir = path(x)
    mkpath(outdir) 

    embed_df = _embed_dict_to_df(x.embedding)

    metadata = Dict(
        :representations => _representations_to_dict(SpeciesRepresentations[x])
    )
    _write_json(joinpath(outdir, "metadata.json"), metadata)
    CSV.write(joinpath(outdir, "representation.csv"), embed_df) 
end

function _load_species_representation(path)
    # todo: gotta change the path save in species rep is relative to projectdir,
    # and isn't absolute

    # okay this is use both by `batch_fit` and by `load_species_representation`
    # (which itself is used by `load_classification_fit`). 
    
    # It might be better if there were two different functions for this. 
    
    # It also must able to deal with multiple cases:
    #. 1. Running `batch_fit` locally with representations generated either
    #locally or on the cluster 
    #  2. Loading batch fits where the batch fit is run on either local/cluster
    #     and the representation is generated on either. 

    is_relative = path[1] == "." 




    #relative_path = is_relative ? "."*string(split(path, "ColoradoBumblebees")[2]) : split(path, "ColoradoBumblebees")
    @info "Top"
    @info split(path, "ColoraodBumblebees")
    relative_path = "."*string(split(path, "ColoradoBumblebees")[2])
    full_path = is_relative ? joinpath(artifactdir()) : joinpath(projectdir(), relative_path) 

    @info is_relative 
    #@info string(split(path, "ColoradoBumblebees")[2]) 
    @info full_path, relative_path, path
    @info "Bottom\n"
    metadata = _read_json(joinpath(full_path, "metadata.json"))
    embed_model = _reconstruct_representation(metadata["representations"])[1]
    embed_dict = _embed_df_to_embed_dict(CSV.read(joinpath(full_path, "representation.csv"), DataFrame)) 
    SpeciesRepresentations(embed_model, embed_dict, relative_path)
end

function _representations_to_dict(reps::Vector{S}) where S<:SpeciesRepresentations
    d = Dict()
    for r in reps
        dict = typeof(r.embed_model) <: RecurrentAutoencoder ? struct2dict(_deserialize(r)) : struct2dict(r.embed_model)
        merge!(d, Dict(string(typeof(r.embed_model)) => dict))
    end
    d
end

function _deserialize(rep::SpeciesRepresentations{RecurrentAutoencoder{V}}) where V
    # convert rep.embed_model to a dict, replace the unserialiables with Symbols
    # of the name, voila 
    d = Dict([f=>getfield(rep.embed_model, f) for f in fieldnames(typeof(rep.embed_model))])

    d[:opt] = Symbol("ADAM_η=$(d[:opt].eta)")
    d[:unit] = Symbol(string(d[:unit]))
    

    return RecurrentAutoencoder{V}(;d...)
end 

function _representations_to_dict(rep::SpeciesRepresentations{T}) where T 
    return T <: RecurrentAutoencoder ? struct2dict(_deserialize(rep).embed_model) : struct2dict(rep.embed_model)
end

function _embed_df_to_embed_dict(embed_df)
    data = load_data()
    dict = Dict{Species,Vector}()
    for r in eachrow(embed_df)
        sp = string(r.species)
        spobj = contains(sp, "Bombus") ? bee(data,sp) : plant(data, sp)
        merge!(dict, Dict{Species,Vector}(spobj => Array(r[1:end-1])))
    end
    dict
end

function _embed_dict_to_df(embed_dict)
    num_species = length(keys(embed_dict))
    per_species_dim = length(embed_dict[keys(embed_dict) |> first])
    cols = []
    push!(cols, [Symbol(fi) => zeros(num_species) for fi in 1:per_species_dim])

    
    df = DataFrame(cols... )
    df.species =  ["" for _ in 1:num_species]
    r = 1
    for (k,v) in embed_dict
        df[r,end] = k.name
        df[r,1:end-1] .= v
        r += 1
    end
    df
end

_reconstruct_representation(::Type{S}, θ) where S<:ColoradoBumblebees.EmbeddingType = nothing

function _reconstruct_representation(representation_metadata)
    # TODO it will be easier to combine this with the save artifacts 
    # SpeciesRepresentations must contain a field that includes its saved path,
    # initialized to nothing initially. 
    # any time `batch_fit` gets called, if the representations aren't saved,
    # make sure they get saved. 

    reps = _get_representation_obj.(keys(representation_metadata))
    
    θs = collect(values(representation_metadata))

    # add dispatch to fix any model specific issues, where fields have to be
    # replaced with types. e.g. for recurrentautoencodder, the opt field needs
    # to be replaced with the object, to just the symbol

    for (i,θ) in enumerate(θs)
        _reconstruct_representation(reps[i], θ)
    end

    reconstructed_reps = [r(; Dict([Symbol(k) => v for (k,v) in θs[i]])...)     for (i,r) in enumerate(reps)]
end


_get_representation_obj(model_string) = Dict(
    "MetawebSVD" => MetawebSVD, 
    "MetawebNode2Vec" => MetawebNode2Vec,
    "LFSVD" => LFSVD, 
    "SimulatedTraits" => SimulatedTraits, 
    "PhylogeneticNode2Vec"=>PhylogeneticNode2Vec,
    "KMeansEnvironmentEmbedding"=>KMeansEnvironmentEmbedding,
    "Pooled"=>Pooled,
    "DenseAutoencoder"=>DenseAutoencoder,
    "RecurrentAutoencoder{Standard}"=>RecurrentAutoencoder{Standard},
    "RecurrentAutoencoder{Variational}"=>RecurrentAutoencoder{Variational},
    "GraphAutoencoder{Standard}" => GraphAutoencoder{Standard},
    "GraphAutoencoder{Variational}" => GraphAutoencoder{Variational}
    )[model_string]


