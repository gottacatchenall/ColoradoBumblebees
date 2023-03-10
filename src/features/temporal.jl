abstract type AutoencoderType end
struct Standard <: AutoencoderType end
struct Variational <: AutoencoderType end

Base.@kwdef struct Autoencoder{T} <: Temporal
    unit = RNN               # node unit for first layer
    encoder_dims = [TEMPORAL_INPUT_DIM, 64, 32]
    decoder_dims = [32, 64, TEMPORAL_INPUT_DIM]
    opt = ADAM(1e-2)
    dropout = 0.1           # dropout
    n_epochs = 1_000        # num epoachs
    train_proportion = 0.8  # 
    batch_size = 64         # batch_size
end
outdim(ae::Autoencoder) = ae.encoder_dims[end]
outdim(ae::Autoencoder, ::Union{Type{Bee},Type{Plant}}) = outdim(ae)

function getfeatures(ae::Autoencoder{Standard}, data)
    enc, dec = _fitmodel(ae, _test_train_split(data, ae.train_proportion)...)
    phen = phenology(data)
    dict = Dict()

    for (k, v) in phen
        merge!(dict, Dict(k => enc(Float32.(v))))
    end
    return dict
end
function getfeatures(ae::Autoencoder{Variational}, data)
    enc_μ, enc_logvar, dec = _fitmodel(ae, _test_train_split(data, ae.train_proportion)...)
    phen = phenology(data)
    dict = Dict()

    for (k, v) in phen
        _, μ, _ = _reparam_trick(enc_μ, enc_logvar, dec, Float32.(v))
        merge!(dict, Dict(k => μ))
    end
    return dict
end

function _fitmodel(ae::Autoencoder{Standard}, train, test)
    enc, dec = _makemodel(ae)
    model = Chain(enc, dec)
    loss(x) = begin
        #Flux.mse(model(x), x)
        mean([Flux.mse(model(t), t) for t in x])
    end

    trainloss, testloss = _train_model(
        model, loss, train, test, ae.opt, ae.n_epochs, ae.batch_size
    )
    return enc, dec
end

function _fitmodel(ae::Autoencoder{Variational}, train, test)
    enc_μ, enc_logvar, dec = _makemodel(ae)
    # model = Chain(enc_μ, enc_logvar, dec)
    function loss(x)
        reconst_loss = 0.0
        kl_div_sum = 0.0
        for t in x
            x̂, μ, logvar = _reparam_trick(enc_μ, enc_logvar, dec, t)
            reconst_loss += Flux.mse(x̂, t)
            kl_div = -0.5 * sum(1.0 .+ logvar .- μ .^ 2 .- exp.(logvar))
            kl_div_sum += kl_div
        end
        return reconst_loss + kl_div_sum
    end

    trainloss, testloss = _train_model(
        (enc_μ, enc_logvar, dec), loss, train, test, ae.opt, ae.n_epochs, ae.batch_size
    )
    return enc_μ, enc_logvar, dec
end

function _train_model(model, loss, train, test, opt, n_epochs, batch_size)
    ps = Flux.params(model)
    trainlosses, testlosses = [], []
    progbar = ProgressMeter.Progress(n_epochs)
    @info "Starting training..."
    for epoch in 1:n_epochs
        Flux.train!(loss, ps, [(train)], opt)

        trainloss = loss(train)
        #  testloss = loss(test)
        testloss = NaN
        ProgressMeter.next!(
            progbar; showvalues=[(Symbol("Train Loss"), trainloss), (Symbol("η"), opt.eta)]
        )
        if epoch % 10 == 0
            push!(trainlosses, trainloss)
            #     push!(testlosses, testloss)
        end
    end
    return trainlosses, testlosses
end

function _test_train_split(data, train_proportion=0.8)
    phen = phenology(data)
    timeseries = shuffle(
        vcat([[Float32.(phen[b]) for b in f(data)] for f in [bees, plants]]...)
    )
    i = Int32(floor(train_proportion * length(timeseries)))
    return timeseries[begin:i], timeseries[(i + 1):end]
end

function _reparam_trick(enc_μ, enc_logvar, dec, x)
    μ, logvar = enc_μ(x), enc_logvar(x)
    std = exp.(logvar ./ 2)
    eps = rand(MvNormal(zeros(Float32, length(μ)), vec(std)))
    x̂ = dec(μ .+ eps .* std)
    return x̂, μ, logvar
end

function _makemodel(ae::Autoencoder{Standard})
    enc_layer_sizes, dec_layer_sizes = ae.encoder_dims, ae.decoder_dims

    _dense_unit(in, out) = Dense(in, out)
    _dense_dropout_unit(in, out) = Chain(Dense(in, out), Dropout(ae.dropout))

    f = ae.dropout > 0 ? _dense_dropout_unit : _dense_unit
    enc = Chain(
        ae.unit(enc_layer_sizes[begin] => enc_layer_sizes[begin + 1]),
        [
            f(enc_layer_sizes[i], enc_layer_sizes[i + 1]) for
            i in 2:(length(enc_layer_sizes) - 1)
        ]...,
    )

    dec = Chain(
        [
            Dense(dec_layer_sizes[i], dec_layer_sizes[i + 1]) for
            i in 1:(length(dec_layer_sizes) - 1)
        ]...,
    )
    return enc, dec
end

function _makemodel(ae::Autoencoder{Variational})
    enc_layer_sizes, dec_layer_sizes = ae.encoder_dims, ae.decoder_dims
    num_enc_layers, num_dec_layers = length.([enc_layer_sizes, dec_layer_sizes])

    _dense_unit(in, out) = Dense(in, out)
    _dense_dropout_unit(in, out) = Chain(Dense(in, out), Dropout(ae.dropout))

    f = ae.dropout > 0 ? _dense_dropout_unit : _dense_unit

    enc_features = Chain(
        ae.unit(enc_layer_sizes[begin] => enc_layer_sizes[begin + 1]),
        [f(enc_layer_sizes[i], enc_layer_sizes[i + 1]) for i in 2:(num_enc_layers - 2)]...,
    )

    encoder_μ = Chain(enc_features, Dense(enc_layer_sizes[end - 1], enc_layer_sizes[end]))
    encoder_logvar = Chain(
        enc_features, Dense(enc_layer_sizes[end - 1], enc_layer_sizes[end])
    )
    dec = Chain(
        [
            Dense(dec_layer_sizes[i], dec_layer_sizes[i + 1]) for i in 1:(num_dec_layers - 1)
        ]...,
    )

    return encoder_μ, encoder_logvar, dec
end
