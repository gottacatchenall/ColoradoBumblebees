using DrWatson
@quickactivate "ColoradoBumblebees"

using DataFrames, CSV
using Random
using Flux
using Flux: DataLoader
using ProgressMeter

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

phen = phenology(data)

label_df = label_dataframe(data)

ts = [Float32.(v) for (k, v) in phen]

enc = Chain(LSTM(TEMPORAL_INPUT_DIM, 64), Dense(64, 32))

dec = Chain(Dense(2 * 32, 32), Dense(32, 8), Dense(8, 1))
model = Chain(x -> (enc(x[1]), enc(x[2])), vcat, dec)

model((ts[1], ts[2]))

features, labels = [], []

for r in eachrow(label_df)
    push!(features, (Float32.(phen[r.bee]), Float32.(phen[r.plant])))
    push!(labels, Float32(r.interaction))
end

loss(x, y) = Flux.mse(first.(model.(x)), y)

features, labels

cut = Int32(floor(0.8length(labels)))
Ishuffled = shuffle(1:length(labels))
Itrain, Itest = Ishuffled[begin:cut], Ishuffled[(cut + 1):end]

batch_size = 256
loader = DataLoader((features[Itrain], labels[Itrain]); batchsize=batch_size)

opt = ADAM()
ps = Flux.params(model)
nepochs = 100

progbar = ProgressMeter.Progress(nepochs)
for e in 1:nepochs
    batchloss = 0.0
    for (x, y) in loader
        Flux.train!(loss, ps, [(x, y)], opt)
        batchloss += loss(x, y)
    end
    ProgressMeter.next!(
        progbar; showvalues=[(Symbol("Batch Loss"), batchloss / batch_size)]
    )
end
