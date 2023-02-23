struct GaussianTemporalEmbedding <: Temporal

end

# this is the extent of min(doy) to max(doy), therefore always the sequence length
const INPUT_DIM = 147


# IDEA:  AE with stacked input that gives matrix embedding for all species at same time 
        # One encoder for everything, one decoder for each species

# Noised and unnoised versions of each 

@Base.kwdef struct AETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::AETemporalEmbedding) = ae.encoder_dims[end]


function getfeatures(ae::AETemporalEmbedding, abundance_timeseries)
    
end

function makemodel(ae::AETemporalEmbedding)
    enc_dims = ae.encoder_dims    
    dec_dims = ae.decoder_dims

    enc = Chain(
        [Dense(enc_dims[i], enc_dims[i+1]) for i in 1:(length(enc_dims)-1)]...
    )
    dec = Chain(
        [Dense(dec_dims[i], dec_dims[i+1]) for i in 1:(length(dec_dims)-1)]...
    )
    enc, dec, Chain(enc,dec) 
end



# ===============


@Base.kwdef struct RNNAETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::RNNAETemporalEmbedding) = ae.encoder_dims[end]

function makemodel(ae::RNNAETemporalEmbedding)
    enc_dims = ae.encoder_dims    
    dec_dims = ae.decoder_dims

    enc = Chain(
        RNN(enc_dims[1] => enc_dims[2]),
        [Dense(enc_dims[i], enc_dims[i+1]) for i in 2:(length(enc_dims)-1)]...
    )
    dec = Chain(
        [Dense(dec_dims[i], dec_dims[i+1]) for i in 1:(length(dec_dims)-1)]...,
    )
    enc, dec, Chain(enc,dec) 
end

# ====================================

@Base.kwdef struct LSTMAETemporalEmbedding <: Temporal
    encoder_dims = [INPUT_DIM, 32, 16]
    decoder_dims = [16, 32, INPUT_DIM]
end
outdim(ae::LSTMAETemporalEmbedding) = ae.encoder_dims[end]

function makemodel(ae::LSTMAETemporalEmbedding)
    enc_dims = ae.encoder_dims    
    dec_dims = ae.decoder_dims

    enc = Chain(
        LSTM(enc_dims[1] => enc_dims[2]),
        [Dense(enc_dims[i], enc_dims[i+1]) for i in 2:(length(enc_dims)-1)]...
    )
    dec = Chain(
        [Dense(dec_dims[i], dec_dims[i+1]) for i in 1:(length(dec_dims)-1)]...,
    )
    enc, dec, Chain(enc,dec) 
end

# ====================================

function fitmodel!(model, timeseries; lr=1e-3, n_epochs=1_000)
    loss(x,y) = mean([Flux.mse(model(t), t) for t in timeseries])


    # X = to_features(timeseries)
    #Itest = sample(1:size(X)[2], Int32(floor(0.2*size(X,2))), replace=false)
    #Itrain = filter(x->x ∉ Itest,  1:size(X,2))
    #data_test = (X[:,Itest],X[:,Itest])
    #data_train = (X[:,Itrain],X[:, Itrain])

    data_test, data_train = (timeseries,timeseries), (timeseries, timeseries)

    data_train = [(ts,ts) for ts in timeseries]

    ps = Flux.params(model)
    
    trainlosses = []
    progbar = ProgressMeter.Progress(n_epochs)
    
    opt = ADAM(lr) 
    
    for epoch in 1:n_epochs

        for (i,ts) in enumerate(timeseries)
            Flux.train!(loss, ps, [data_train[i]], opt)            
        end 
    
        trainloss = mean([loss(dt...) for dt in data_train])
        ProgressMeter.next!(
                progbar;
                showvalues = [
                    (Symbol("Loss"), trainloss),
        
                ],
            )
        if epoch % 10 == 0
            push!(trainlosses, trainloss)
        end
    end
end


# tests

phen = phenology(data)

bees_timeseries = [Float32.(phen[b]) for b in bees(data)]

embeds = []

enc, dec, mod = makemodel(RNNAETemporalEmbedding())
fitmodel!(mod, bees_timeseries)



enc, dec, mod = makemodel(LSTMAETemporalEmbedding())
fitmodel!(mod, bees_timeseries)


enc, dec, mod = makemodel(AETemporalEmbedding())
fitmodel!(mod, bees_timeseries)

enc.(bees_timeseries)

f = Figure()
ax1 = Axis(f[1,1], title="True")
ax2 = Axis(f[2,1], title="Reconstructed")


scatterlines!(ax1, 
    1:INPUT_DIM,
     [sum([bees_timeseries[i][t] for i in 1:length(bees_timeseries)]) for t in 1:INPUT_DIM]
)
scatterlines!(ax1, 
    1:INPUT_DIM,
     [sum([mod.(bees_timeseries[i])[t] for i in 1:length(bees_timeseries)]) for t in 1:INPUT_DIM]
)
current_figure()


f = Figure()
ax1 = Axis(f[1,1])
for v in enc.(bees_timeseries)
    scatter!(ax1, [v[1]], [v[2]])
end


# Functions to construct inputs:
# There are diffrent ways to do this:
#       1. pooled 
#       2. site stratified
#       3. year stratified
#       4. year & site strat


function phenology(data)    
    species = vcat(bees(data)..., plants(data)...)
    firstdoy, lastdoy = extrema([dayofyear(i.time) for i in interactions(data)])

    dict = Dict()
    for sp in species
        abundances = zeros(Int32, lastdoy-firstdoy+1)
        ints = interactions(data, sp)
        for i in ints
            abundances[dayofyear(i.time) - firstdoy + 1] += 1
        end
        merge!(dict, Dict(sp=>abundances))
    end
    dict
end