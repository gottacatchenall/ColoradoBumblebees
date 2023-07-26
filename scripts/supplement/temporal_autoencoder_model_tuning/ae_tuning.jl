using CUDA
using Flux

using DataFrames
using CSV
using ProgressMeter

function load_data()
    df = CSV.read("./pheno.csv", DataFrame)

    timeseries = []

    species = unique(df.species)
    for s in species
        this_species = filter(r->r.species == s, df)
        I = sortperm(this_species.time)
        push!(timeseries, this_species.abundance[I])
    end

    species, timeseries
end 

function get_data_loader()
    _, timeseries = load_data()
    train_loader = Flux.DataLoader((timeseries, timeseries), batchsize=length(timeseries), shuffle=true)
end

function _makemodel(unit, rnn_dims, encoder_dims, decoder_dims)
    rnn = _make_rnn_chain(unit, rnn_dims)
    enc = _make_dense_chain(encoder_dims)
    dec = _make_dense_chain(decoder_dims)
    return rnn, enc, dec
end

function _make_rnn_chain(unit, dims_per_layer)
    Chain([unit(dims_per_layer[i]=>dims_per_layer[i+1]) for i in 1:(length(dims_per_layer)-1)]...)
end

function _make_dense_chain(dims_per_layer)
    Chain([Dense(dims_per_layer[i], dims_per_layer[i + 1]) for i in 1:(length(dims_per_layer)-1)]...)
end

function loss(rnn, enc, dec, x)
    reconst_loss = 0f0
    for t in x
        Flux.reset!(rnn) # Reset hidden state
        hidden = vcat([rnn([ti]) for ti in t]...)
        reconst_loss += Flux.mse(dec(enc(vec(hidden))), t)
    end
    return reconst_loss
end

function train(unit, rnn_dims, encoder_dims, decoder_dims; η=0.01, n_epochs=1000, cuda=false)
    loader = get_data_loader()

    rnn, enc, dec = _makemodel(unit, rnn_dims, encoder_dims, decoder_dims)

    if cuda 
        rnn |> gpu 
        enc |> gpu
        dec |> gpu 
        loader |> gpu 
    end

    opt = ADAM(η)
    ps = Flux.params(rnn, enc, dec)
    progbar = ProgressMeter.Progress(2*n_epochs)
    losses = []
    for _ in 1:n_epochs
        for (x,_) in loader
            Flux.train!(x->loss(rnn,enc,dec,x), ps, [(x)], opt)
            trainloss = loss(rnn, enc, dec, x)
            push!(losses, trainloss)
            ProgressMeter.next!(
                progbar; showvalues=[(Symbol("Train Loss"), trainloss)]
            )
        end
    end

    rnn, enc
end

function write_embeddings(rnn, enc, filename; cuda=false)
    species, timeseries = load_data()
    
    embeds = []

    if cuda 
        timeseries |> gpu
    end 
    for t in timeseries
        Flux.reset!(rnn) # Reset hidden state
        hidden = vcat([rnn([ti]) for ti in t]...)
        out = enc(hidden) |> cpu
        push!(embeds, out)
    end
    
    mat = hcat(embeds...)
    
    df = DataFrame(Matrix(mat'), :auto)

    df.species = species

    CSV.write("./$filename.csv", df)
end

function main()
    cuda = true
    rnn, enc = train(LSTM, [1, 8, 1], [147,  16], [16, 147]; cuda=cuda)
    write_embeddings(rnn,enc,filename; cuda=cuda)
end

