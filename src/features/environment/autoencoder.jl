@Base.kwdef struct EnvironmentAutoencoder <: Environment
    species_encoder_dims = [16]
    species_decoder_dims = [128, 256]
    shared_decoder_dims = [16, 32, 64, 128]
    n_epochs = 1_000
    optimizer = ADAM(1e-3)
    full = false
    obs_size = 15000
end
outdim(ae::EnvironmentAutoencoder) = ae.species_encoder_dims[end]
outdim(ae::EnvironmentAutoencoder, ::Union{Type{Bee},Type{Plant}}) = outdim(ae)

function getfeatures(ae::EnvironmentAutoencoder, data)
    encoders, xs = _fitmodel(ae, data)
    dict = Dict()
    for (k,m) in encoders
        merge!(dict, Dict(k=>m[:encoder](xs[k])))
    end
    dict
end

function _fitmodel(ae::EnvironmentAutoencoder, data)
    env = CSV.read(datadir("public","environment", "whitened.csv"), DataFrame)
    
    species = vcat(bees(data), plants(data))
    species_dfs = _get_species_dfs(ae, species, env)

    xs = Dict([s=>Float32.(vec(Matrix(species_dfs[s][!,2:end]))) for s in species])
    
    species_specific, shared = _makemodel(ae, species, species_dfs)

    models = [Chain(species_specific[s][:encoder], shared, species_specific[s][:decoder]) for s in species]
    losses = [(x,y) -> Flux.mse(models[i](x), y) for i in eachindex(models)]

    opt = ae.optimizer
    ps = Flux.params(models)
    progbar = ProgressMeter.Progress(ae.n_epochs)
    for _ in 1:ae.n_epochs
        s = 0.
        for (i,sp) in enumerate(species)
            x = xs[sp]
            if length(x) > 30
                Flux.train!(losses[i], ps, [(x,x)], opt)         
                trainloss = Flux.mse(models[i](x), x)
                s += trainloss
            end 
        end
        ProgressMeter.next!(progbar; showvalues = [(Symbol("Train Loss"), s)])
    end    

    species_specific, xs
end

function _makemodel(ae::EnvironmentAutoencoder, species, species_dfs)
    _species_specific_models(ae, species, species_dfs), _shared_decoder(ae)
end

function _species_specific_models(ae::EnvironmentAutoencoder, species, species_dfs)
    species_models = Dict()

    for s in species
        enc, dec = _species_specific_encoder_decoder(ae, species_dfs[s])
        merge!(species_models, Dict(s=>Dict(:encoder=>enc, :decoder=>dec)))
    end
    
    species_models
end

_filter_species(env_df, species) = filter(species, env_df)

function _shrink(ae, env)
    total_obs = ae.obs_size
    beedf = filter(x-> split(x.species, " ")[1] == "Bombus", env)

    rowI = shuffle(1:nrow(beedf))[1:(Int32(floor(total_obs/2)))]
    beedf = beedf[rowI,:]
    plantdf = filter(x-> split(x.species, " ")[1] != "Bombus", env)
    vcat(plantdf[shuffle(1:nrow(plantdf))[1:nrow(beedf)],:], beedf)
end

function  _get_species_dfs(ae, species, env_df)
    df = env_df
    if !ae.full
        @info "You are running a smaller environmental embedding with $(ae.obs_size) observations for both plants/bees.\nThis is changed via the kwarg full=true"
        df = _shrink(ae, env_df)
    end

    species_dfs = Dict()
    for s in species
        thisdf = filter(x-> x.species == s.name, df)
        merge!(species_dfs, Dict(s=>thisdf))
    end
    species_dfs
end

function _species_specific_encoder_decoder(ae, this_species_env_df)
    concat = vec(Matrix(this_species_env_df[!,2:end]))

    enc = Chain(
        Dense(length(concat), ae.species_encoder_dims[begin]),
        [Dense(ae.species_encoder_dims[t], ae.species_encoder_dims[t+1]) for t in 1:length(ae.species_encoder_dims)-1]...,
    )
        
    dec = Chain(
        [Dense(ae.species_decoder_dims[t], ae.species_decoder_dims[t+1]) for t in 1:length(ae.species_decoder_dims)-1]...,
        Dense(ae.species_decoder_dims[end], length(concat)),
        )
    enc,dec
end

function _shared_decoder(ae::EnvironmentAutoencoder)
    shared_dec = Chain(
        [Dense(ae.shared_decoder_dims[t], ae.shared_decoder_dims[t+1]) for t in 1:length(ae.shared_decoder_dims)-1]...,
    )
end