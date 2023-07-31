Base.@kwdef struct RecurrentAutoencoder <: Temporal
    unit = RNN               # node unit for first layer
    rnn_dims = [1, 8, 1]
    encoder_dims = [TEMPORAL_INPUT_DIM, 8]
    decoder_dims = [8, TEMPORAL_INPUT_DIM]
    opt = ADAM(1e-2)
    n_epochs = 500        # num epoachs
end
outdim(ae::RecurrentAutoencoder) = ae.encoder_dims[end]
outdim(ae::RecurrentAutoencoder, ::Union{Type{Bee},Type{Plant}}) = outdim(ae)

function getfeatures(ae::RecurrentAutoencoder, data)
    phen = load_phenology(data)
    rnn, enc, _ = _fitmodel(ae, phen)

    dict = Dict()

    for (k, v) in phen
        Flux.reset!(rnn) # Reset hidden state
        hidden = vcat([rnn([ti]) for ti in v]...)
        emb = enc(hidden)
        merge!(dict, Dict(k => emb))
    end
    return dict
end

function _fitmodel(ae::RecurrentAutoencoder, phenologies)
    rnn, enc, dec = _makemodel(ae)
    _train_model!(ae, rnn, enc, dec, phenologies)
    return rnn, enc, dec
end


function _train_model!(ae::RecurrentAutoencoder, rnn, enc, dec, phenologies)

    function loss(x)
        reconst_loss = 0f0
        for t in values(phenologies)
            Flux.reset!(rnn) # Reset hidden state
            hidden = vcat([rnn([ti]) for ti in t]...)
            #concat_hidden = vec(hcat(hidden...))
            reconst_loss += Flux.mse(dec(enc(hidden)), t)
        end
        return reconst_loss
    end

    opt, n_epochs = ae.opt, ae.n_epochs
    ps = Flux.params(rnn, enc, dec)
    progbar = ProgressMeter.Progress(n_epochs)
    @info "Starting training..."
    for epoch in 1:n_epochs
        #=∇ = gradient(ps) do 
            # Compute MSE loss on rest of sequence
            ls = 0f0
            for t in values(phenologies)
                Flux.reset!(rnn) # Reset hidden state
                hidden = vcat([rnn([ti]) for ti in t]...)
                #concat_hidden = vec(hcat(hidden...))
                ls += Flux.mse(dec(enc(hidden)), t)
            end
               
            ls 
        end
        ProgressMeter.next!(progbar)
        Flux.update!(opt, ps, ∇)=#

        Flux.train!(loss, ps, [(phenologies)], opt)
        trainloss = loss(phenologies)
        ProgressMeter.next!(
            progbar; showvalues=[(Symbol("Train Loss"), trainloss)]
        )

    end
end

function _makemodel(ae::RecurrentAutoencoder)
    
    rnn = _make_rnn_chain(ae.unit, ae.rnn_dims)
    enc = _make_dense_chain(ae.encoder_dims)
    dec = _make_dense_chain(ae.decoder_dims)

    return rnn, enc, dec
end

function _make_rnn_chain(unit, dims_per_layer)
    Chain([unit(dims_per_layer[i]=>dims_per_layer[i+1]) for i in 1:(length(dims_per_layer)-1)]...)
end

function _make_dense_chain(dims_per_layer)
    Chain([Dense(dims_per_layer[i], dims_per_layer[i + 1]) for i in 1:(length(dims_per_layer)-1)]...)
end