
# IDEAS:
#  -   preprocess signal with local smoothing (large fluctuations are likely
#      measurement error and also could create spikier loss surface) 

# - Each Temporal model is a struct with follwing params
#       - enc & dec dims
#       - node unit (dense/nothing, RNN, LSTM, GRU)
#       - variational
#       - noisy
#
# so only need one struct for AEs and rest is handled later

# this is the extent of min(doy) to max(doy), therefore always the sequence length
const INPUT_DIM = 147

Base.@kwdef struct AETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::AETemporalEmbedding) = ae.encoder_dims[end]
outdim(ae::AETemporalEmbedding, ::Union{Type{Bee},Type{Plant}}) = outdim(ae)

function makemodel(ae::AETemporalEmbedding)
    enc_dims = ae.encoder_dims
    dec_dims = ae.decoder_dims

    enc = Chain([Dense(enc_dims[i], enc_dims[i + 1]) for i in 1:(length(enc_dims) - 1)]...)
    dec = Chain([Dense(dec_dims[i], dec_dims[i + 1]) for i in 1:(length(dec_dims) - 1)]...)
    return enc, dec, Chain(enc, dec)
end

# ======================

Base.@kwdef struct VAETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 128, 64, 32]
    decoder_dims = [32, 64, INPUT_DIM]
end
outdim(vae::VAETemporalEmbedding) = vae.encoder_dims[end]
function makemodel(vae::VAETemporalEmbedding)
    enc_dims, dec_dims = vae.encoder_dims, vae.decoder_dims

    enc_features = Chain(
        LSTM(enc_dims[begin] => enc_dims[2]),
        [Dense(enc_dims[i], enc_dims[i + 1]) for i in 2:(length(enc_dims) - 2)]...,
    )
    encoder_μ = Chain(enc_features, Dense(enc_dims[end - 1], enc_dims[end]))
    encoder_logvar = Chain(enc_features, Dense(enc_dims[end - 1], enc_dims[end]))
    dec = Chain([Dense(dec_dims[i], dec_dims[i + 1]) for i in 1:(length(dec_dims) - 1)]...)

    return encoder_μ, encoder_logvar, dec
end

# ===============

Base.@kwdef struct RNNAETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::RNNAETemporalEmbedding) = ae.encoder_dims[end]

function makemodel(ae::RNNAETemporalEmbedding)
    enc_dims = ae.encoder_dims
    dec_dims = ae.decoder_dims
    enc = Chain(
        RNN(enc_dims[1] => enc_dims[2]),
        [Dense(enc_dims[i], enc_dims[i + 1]) for i in 2:(length(enc_dims) - 1)]...,
    )
    dec = Chain([Dense(dec_dims[i], dec_dims[i + 1]) for i in 1:(length(dec_dims) - 1)]...)

    return enc, dec, Chain(enc, dec)
end

# ====================================

Base.@kwdef struct LSTMAETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::LSTMAETemporalEmbedding) = ae.encoder_dims[end]

function makemodel(ae::LSTMAETemporalEmbedding)
    enc_dims = ae.encoder_dims
    dec_dims = ae.decoder_dims

    enc = Chain(
        LSTM(enc_dims[1] => enc_dims[2]),
        [Dense(enc_dims[i], enc_dims[i + 1]) for i in 2:(length(enc_dims) - 1)]...,
    )
    dec = Chain([Dense(dec_dims[i], dec_dims[i + 1]) for i in 1:(length(dec_dims) - 1)]...)
    return enc, dec, Chain(enc, dec)
end

# ====================================

function fitmodel!(model, timeseries; lr=1e-3, n_epochs=1_000, noise=false, σ=0.03)
    data_train = if noise
        [(ts .+= rand(Normal(0, σ)), ts) for ts in timeseries]
    else
        [(ts, ts) for ts in timeseries]
    end

    loss(x, y) = mean([Flux.mse(model(t), t) for t in timeseries])
    ps = Flux.params(model)

    trainlosses = []
    progbar = ProgressMeter.Progress(n_epochs)

    opt = ADAM(lr)

    for epoch in 1:n_epochs
        for (i, ts) in enumerate(timeseries)
            in, out = data_train[i]
            Flux.train!(loss, ps, [(in, out)], opt)
        end

        trainloss = mean([loss(dt...) for dt in data_train])
        ProgressMeter.next!(progbar; showvalues=[(Symbol("Loss"), trainloss)])
        if epoch % 10 == 0
            push!(trainlosses, trainloss)
        end
    end
end

using Flux
using Flux: logitbinarycrossentropy
using ProgressMeter: Progress, next!
using Zygote

# tests
phen = phenology(data)

bees_timeseries = [Float32.(phen[b]) for b in bees(data)]

embeds = []

enc, dec, mod = makemodel(RNNAETemporalEmbedding([INPUT_DIM, 64, 32], [32, 64, INPUT_DIM]))
fitmodel!(mod, bees_timeseries; lr=1e-3, noise=false)

enc, dec, mod = makemodel(
    LSTMAETemporalEmbedding([INPUT_DIM, 256, 128, 64, 32], [32, 64, 128, INPUT_DIM])
)
fitmodel!(mod, bees_timeseries; lr=1e-4)

enc, dec, mod = makemodel(AETemporalEmbedding())
fitmodel!(mod, bees_timeseries)

# ====================================
# 
#   VAE 
# 
# ============================================

function reparam(model, x)
    enc_μ, enc_logvar, dec = model
    μ, logvar = enc_μ(x), enc_logvar(x)
    std = exp.(logvar ./ 2)
    eps = rand(MvNormal(zeros(Float32, length(μ)), vec(std)))
    x̂ = dec(μ .+ eps .* std)
    return x̂, μ, logvar
end

function fitvae!(model, timeseries; lr=1e-3, n_epochs=1_000)
    data_train = [(ts, ts) for ts in timeseries]

    function loss(x, y)
        x̂, μ, logvar = reparam(model, x)
        reconst_loss = sum([Flux.mse(x̂, t) for t in timeseries])

        kl_div = -0.5 * sum(1.0 .+ logvar .- μ .^ 2 .- exp.(logvar))
        return reconst_loss + kl_div
    end

    ps = Flux.params(model)

    trainlosses = []
    progbar = ProgressMeter.Progress(n_epochs)

    opt = ADAM(lr)

    for epoch in 1:n_epochs
        for (i, ts) in enumerate(timeseries)
            in, out = data_train[i]

            Flux.train!(loss, ps, [(in, out)], opt)
        end

        trainloss = mean([loss(dt...) for dt in data_train])
        ProgressMeter.next!(progbar; showvalues=[(Symbol("Loss"), trainloss)])
        if epoch % 10 == 0
            push!(trainlosses, trainloss)
        end
    end
end

enc_μ, enc_logvar, decoder = makemodel(VAETemporalEmbedding())

fitvae!((enc_μ, enc_logvar, decoder), bees_timeseries; lr=1e-5, n_epochs=10_000)

xpred = [reparam((enc_μ, enc_logvar, decoder), ts)[1] for ts in bees_timeseries]

f = Figure()
ax1 = Axis(f[1, 1]; title="")
ax2 = Axis(f[1, 2]; title="")

barplot!(
    ax1,
    1:INPUT_DIM,
    [sum([bees_timeseries[i][t] for i in 1:length(bees_timeseries)]) for t in 1:INPUT_DIM],
)
barplot!(
    ax2,
    1:INPUT_DIM,
    [sum([xpred[i][t] for i in 1:length(bees_timeseries)]) for t in 1:INPUT_DIM],
)
current_figure()

scatter(xpred[2] .- bees_timeseries[2])

# ==============================================

f = Figure()
ax1 = Axis(f[1, 1]; title="True")
ax2 = Axis(f[2, 1]; title="Reconstructed")
scatterlines!(
    ax1,
    1:INPUT_DIM,
    [sum([bees_timeseries[i][t] for i in 1:length(bees_timeseries)]) for t in 1:INPUT_DIM],
)
scatterlines!(
    ax2,
    1:INPUT_DIM,
    [
        sum([mod(bees_timeseries[i])[t] for i in 1:length(bees_timeseries)]) for
        t in 1:INPUT_DIM
    ],
)
current_figure()

# =============================================================
f = Figure()
ax1 = Axis(f[1, 1]; title="")
beenum = 2
scatterlines!(ax1, 1:INPUT_DIM, bees_timeseries[beenum])
scatterlines!(ax1, 1:INPUT_DIM, mod(bees_timeseries[beenum]))
current_figure()

f = Figure()
ax1 = Axis(f[1, 1])
for v in enc.(bees_timeseries)
    scatter!(ax1, [v[1]], [v[2]])
end
current_figure()

# Functions to construct inputs:
# There are diffrent ways to do this:
#       1. pooled 
#       2. site stratified
#       3. year stratified
#       4. year & site strat
