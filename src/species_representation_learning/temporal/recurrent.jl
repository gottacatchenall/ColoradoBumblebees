Base.@kwdef struct RecurrentAutoencoder{T<:AutoencoderType} <: Temporal
    unit = RNN               # node unit for first layer
    rnn_dims = [1, 8, 1]
    encoder_dims = [TEMPORAL_INPUT_DIM, 8]
    decoder_dims = [8, TEMPORAL_INPUT_DIM]
    opt = ADAM(1e-2)
    n_epochs = 500       
end
outdim(ae::RecurrentAutoencoder) = ae.encoder_dims[end]
outdim(ae::RecurrentAutoencoder, ::Union{Type{Bee},Type{Plant}}) = outdim(ae)

function DrWatson.savename(ae::RA) where RA<:RecurrentAutoencoder
    η = ae.opt.eta
    
    rnn_dims = ae.rnn_dims
    dec_dims = ae.decoder_dims
    unit = string(ae.unit)
    @info ae
    @info enc_dims, dec_dims

    "unit_$(unit)_learningrate_$(η)_rnn_$(rnn_dims)_dec_$(dec_dims)_nepochs_$(ae.n_epochs)"
end

function _embed(data::BeeData, ae::RecurrentAutoencoder{Standard})
    phen = load_phenology(data)
    rnn, enc, _ = _fitmodel(ae, phen)

    dict = Dict()

    for (k, v) in phen
        Flux.reset!(rnn) # Reset hidden state
        hidden = [rnn([ti])[1] for ti in v]
        emb = enc(hidden)
        merge!(dict, Dict(k => emb))
    end
    return dict
end

function _embed(data::BeeData, ae::RecurrentAutoencoder{Variational})
    phen = load_phenology(data)
    rnn, first_enc, enc_μ, enc_logvar, dec = _fitmodel(ae, phen)

    dict = Dict()

    for (k, v) in phen
        Flux.reset!(rnn) # Reset hidden state
        hidden = [rnn([ti])[1] for ti in v]
        emb = enc_μ(first_enc(hidden))
        merge!(dict, Dict(k => emb))
    end
    return dict
end

function _get_data_loader(phenlogies)
    timeseries = collect(values(phenlogies))    
    return Flux.DataLoader((timeseries, timeseries), batchsize=length(timeseries), shuffle=true)
end

convert_to_gpu!(x) = x |> gpu 

function _fitmodel(ae::RecurrentAutoencoder{Standard}, phenologies)
    loader = _get_data_loader(phenologies)
    rnn, enc, dec = _makemodel(ae)
    
    ColoradoBumblebees.GPU_AVAILABLE && convert_to_gpu!.([rnn, enc, dec, loader])

    _train_model!(ae, rnn, enc, dec, loader)

    return rnn, enc, dec
end

function _fitmodel(ae::RecurrentAutoencoder{Variational}, phenologies)
    loader = _get_data_loader(phenologies)
    rnn, first_enc, enc_μ, enc_logvar, dec = _makemodel(ae)

    ColoradoBumblebees.GPU_AVAILABLE && convert_to_gpu!.([rnn, first_enc, enc_μ, enc_logvar, dec, loader])

    _train_model!(ae, rnn, first_enc, enc_μ, enc_logvar, dec, loader)

    return rnn, first_enc, enc_μ, enc_logvar, dec
end


function _train_model!(ae::RecurrentAutoencoder{Variational}, rnn, first_enc, enc_μ, enc_logvar, dec, loader)
    function _reparam_trick(first_enc, enc_μ, enc_logvar, dec, x)
        μ, logvar = enc_μ(first_enc(x)), enc_logvar(first_enc(x))
        std = exp.(logvar ./ 2)
        eps = rand(MvNormal(zeros(Float32, length(μ)), vec(std)))
        x̂ = dec(μ .+ eps .* std)
        return x̂, μ, logvar 
    end
    function variational_loss(rnn, first_enc, enc_μ, enc_logvar,dec,x)
        reconst_loss = 0f0
        kl_div_sum = 0f0
        for t in values(x)
            Flux.reset!(rnn) # Reset hidden state
            hidden = [rnn([ti])[1] for ti in t]
            x̂, μ, logvar = _reparam_trick(first_enc, enc_μ, enc_logvar, dec, hidden)
            reconst_loss += Flux.mse(x̂, t)
            kl_div = -0.5 * sum(1.0 .+ logvar .- μ .^ 2 .- exp.(logvar))
            kl_div_sum += kl_div
        end
        return reconst_loss + kl_div_sum
    end

    opt, n_epochs = ae.opt, ae.n_epochs
    ps = Flux.params(rnn,first_enc, dec)
    progbar = ProgressMeter.Progress(n_epochs)
    @info "Starting training..."
    for epoch in 1:n_epochs
        for (x,_) in loader
            Flux.train!(x->variational_loss(rnn, first_enc, enc_μ, enc_logvar,dec, x), ps, [(x)], opt)
            trainloss = variational_loss(rnn, first_enc, enc_μ, enc_logvar,dec, x)
            ProgressMeter.next!(
                progbar; showvalues=[(Symbol("Train Loss"), trainloss)]
            )
        end
    end
end 



function _train_model!(ae::RecurrentAutoencoder{Standard}, rnn, enc, dec, loader)
    function loss(rnn, enc, dec, x)
        reconst_loss = 0f0
        for t in x
            Flux.reset!(rnn) # Reset hidden state
            hidden = vcat([rnn([ti]) for ti in t]...)
            reconst_loss += Flux.mse(dec(enc(vec(hidden))), t)
        end
        return reconst_loss
    end
    

    opt, n_epochs = ae.opt, ae.n_epochs
    ps = Flux.params(rnn, enc, dec)
    progbar = ProgressMeter.Progress(n_epochs)
    @info "Starting training..."
    for epoch in 1:n_epochs
        for (x,_) in loader
            Flux.train!(x->loss(rnn,enc,dec,x), ps, [(x)], opt)
            trainloss = loss(rnn, enc, dec, x)
            ProgressMeter.next!(
                progbar; showvalues=[(Symbol("Train Loss"), trainloss)]
            )
        end
    end
end

function _makemodel(ae::RecurrentAutoencoder{Standard})
    rnn = _make_rnn_chain(ae.unit, ae.rnn_dims)
    enc = _make_dense_chain(ae.encoder_dims)
    dec = _make_dense_chain(ae.decoder_dims)

    return rnn, enc, dec
end


function _makemodel(ae::RecurrentAutoencoder{Variational})
    rnn = _make_rnn_chain(ae.unit, ae.rnn_dims)

    first_enc = _make_dense_chain(ae.encoder_dims[1:end])
    enc_μ = Chain(Dense(ae.encoder_dims[end], ae.encoder_dims[end]))
    enc_logvar = Chain(Dense(ae.encoder_dims[end], ae.encoder_dims[end]))

    dec = _make_dense_chain(ae.decoder_dims)

    return rnn, first_enc, enc_μ, enc_logvar, dec
end

function _make_rnn_chain(unit, dims_per_layer)
    Chain([unit(dims_per_layer[i]=>dims_per_layer[i+1]) for i in 1:(length(dims_per_layer)-1)]...)
end

function _make_dense_chain(dims_per_layer)
    Chain([Dense(dims_per_layer[i], dims_per_layer[i + 1]) for i in 1:(length(dims_per_layer)-1)]...)
end

