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
    metadata = _read_json(joinpath(path, "metadata.json"))
    embed_model = _reconstruct_representation(metadata["representations"])[1]
    embed_dict = _embed_df_to_embed_dict(CSV.read(joinpath(path, "representation.csv"), DataFrame)) 
    SpeciesRepresentations(embed_model, embed_dict, path)
end

function _representations_to_dict(reps::Vector{S}) where S<:SpeciesRepresentations
    d = Dict()
    for r in reps
        merge!(d, Dict(string(typeof(r.embed_model)) => struct2dict(r.embed_model)))
    end
    d
end

function _representations_to_dict(rep::SpeciesRepresentations)
    struct2dict(rep.embed_model)
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

function _reconstruct_representation(representation_metadata)
    # TODO it will be easier to combine this with the save artifacts 
    # SpeciesRepresentations must contain a field that includes its saved path,
    # initialized to nothing initially. 
    # any time `batch_fit` gets called, if the representations aren't saved,
    # make sure they get saved. 

    @info representation_metadata

    reps = _get_representation_obj.(keys(representation_metadata))
    
    θs = collect(values(representation_metadata))

    reconstructed_reps = [r(; Dict([Symbol(k) => v for (k,v) in θs[i]])...)     for (i,r) in enumerate(reps)]
end


_get_representation_obj(model_string) = Dict(
    "MetawebSVD" => MetawebSVD, 
    "LFSVD" => LFSVD, 
    "SimulatedTraits" => SimulatedTraits, 
    "PhylogeneticNode2Vec"=>PhylogeneticNode2Vec,
    "KMeansEnvironmentEmbedding"=>KMeansEnvironmentEmbedding,
    "Pooled"=>Pooled,
    "DenseAutoencoder"=>DenseAutoencoder,
    "RecurrentAutoencoder"=>RecurrentAutoencoder
    )[model_string]


