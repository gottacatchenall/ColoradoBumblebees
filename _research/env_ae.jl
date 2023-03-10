using DrWatson
@quickactivate "ColoradoBumblebees"

using DataFrames, CSV
using Random
using Flux
using Flux: DataLoader
using ProgressMeter
using MLJ
using DataFrames
using Distributions
using ParameterSchedulers
using ParameterSchedulers: Scheduler

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

embeddings = [
    EnvironmentAutoencoder(;
        n_epochs=200,
        species_encoder_dims=[256, 8],
        shared_decoder_dims=[8, 64, 128],
        species_decoder_dims=[128, 256],
    ),
]

df = feature_dataframe(data, embeddings)

DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
BRT = @load EvoTreeClassifier pkg = EvoTrees

rf = RandomForest()
dt = DecisionTree()
brt = BRT()

function single_run(X, y, ens_size=256, batch_size=128)
    Is = shuffle(1:nrow(df))
    cut = Int32(floor(0.8 * nrow(df)))
    Itrain, Itest = Is[1:cut], Is[(cut + 1):end]
    ytest = [x == true for x in y[Itest]]
    ypredict = zeros(length(Itest))
    for i in 1:ens_size
        mach = machine(brt, X, y)
        theserows = balance_sample(y, Itrain, batch_size, 0.5)
        fit!(mach; rows=theserows, verbosity=0)
        pred = predict(mach; rows=Itest)
        ypredict .+= [p.prob_given_ref[2] for p in pred]
    end

    ypredict = ypredict ./ (ens_size)
    return computemeasures(ytest, ypredict)
end

y, X, _, = unpack(df, ==(:interaction), ∉([:bee, :plant]); rng=123)
y = coerce(y, Multiclass{2})

single_run(X, y)

print()

#=
env = CSV.read(datadir("public","environment", "whitened.csv"), DataFrame)

sp = vcat(bees(data), plants(data))

total_obs = 15000
beedf = filter(x-> split(x.species, " ")[1] == "Bombus", env)

rowI = shuffle(1:nrow(beedf))[1:(Int32(floor(total_obs/2)))]
beedf = beedf[rowI,:]
plantdf = filter(x-> split(x.species, " ")[1] != "Bombus", env)

train_df = vcat(plantdf[shuffle(1:nrow(plantdf))[1:nrow(beedf)],:], beedf)

labs = zeros(Float32, length(sp), 2*nrow(beedf))
feat = zeros(Float32, 19, 2*nrow(beedf))

train_df = vcat(plantdf[shuffle(1:nrow(plantdf))[1:nrow(beedf)],:], beedf)

species_dfs = Dict()
cursor = 1
for s in sp
    thisdf = filter(x-> x.species == s.name, train_df)
    merge!(species_dfs, Dict(s=>thisdf))
   #= i = findfirst(isequal(s), sp)
    for r in eachrow(thisdf)
        labs[i, cursor] = 1
        feat[:,cursor] .= Vector(r[2:end])
        cursor += 1
    end =#
end

feat
labs

onehot(s) = [i == findfirst(isequal(s), sp) for (i,_) in enumerate(sp)]

one_hot_dict = Dict([s=>onehot(s) for s in sp])

function make_species_encoder_decoder(s, in_embedding_dim, out_embedding_dim)
    concat = vec(Matrix(species_dfs[s][!,2:end]))

    enc = Chain(
        Dense(length(concat), 64, relu),
        Dense(64, in_embedding_dim))
    dec = Chain(Dense(out_embedding_dim, length(concat)))
    enc,dec
end

embedding_dim = 16
shared_decoder_out_dim = 128

shared_dec = Chain(
    Dense(embedding_dim, 32,relu),
    Dense(32,64,relu),
    Dense(64,shared_decoder_out_dim)
)

# not enough obs of this
sp = filter(s->s.name != "Hypochaeris radicata", sp)

species_enc, species_dec = Dict(), Dict()

for s in sp
    enc, dec = make_species_encoder_decoder(s, embedding_dim, shared_decoder_out_dim)
    merge!(species_enc, Dict(s=>enc))
    merge!(species_dec, Dict(s=>dec))
end

n_epochs = 200

models = [Chain(species_enc[b], shared_dec, species_dec[b]) for b in sp]
ps = Flux.params(models)
xs = [Float32.(vec(Matrix(species_dfs[b][!,2:end]))) for b in sp]
progbar = ProgressMeter.Progress(n_epochs)
for epoch in 1:n_epochs
    s = 0.
    for (i,b) in enumerate(sp)
        x = xs[i]
        Flux.train!(losses[i], ps, [(x,x)], opt)         
        trainloss =Flux.mse(models[i](x), x)
        s += trainloss
    end
    ProgressMeter.next!(progbar; showvalues = [(Symbol("Train Loss"), s)])
end

# cut = Int32(floor(0.8data_size))
# ishuffled = shuffle(1:nrow(train_df))
# Itrain, Itest = ishuffled[1:cut], ishuffled[cut+1:data_size]
#Xtrain, Ytrain = feat[:,Itrain], labs[:,Itrain]
#Xtest, Ytest = feat[:,Itest], labs[:,Itest]
#=
envlayers = 19
numbees = length(sp)

enc = Chain(
    Dense(envlayers, 256, relu),
    Dense(256, 64, relu),
    Dense(64, 8)
)
dec = Chain(
    Dense(8, 32, relu),
    Dense(32, envlayers),
)
model = Chain(enc,dec)

opt = ADAM(1e-2)
ps = Flux.params(model)

loss(x,y) = Flux.mse(model(x),y)

n_epochs = 200
batch_size = 1028

train_loader = DataLoader((Xtrain, Ytrain); batchsize=batch_size)
progbar = ProgressMeter.Progress(n_epochs)
for epoch in 1:n_epochs
    epoch_sum = 0.
    for (x, y) in train_loader 
        Flux.train!(loss, ps, [(x,y)], opt)  
        epoch_sum += loss(x,y)
    end
    ProgressMeter.next!(progbar; showvalues = [(Symbol("Batch Loss"), epoch_sum/batch_size)])
end

enc(Xtest)

=#

length(findall(x->length(x) > 0,[interactions(data,b,p) for b in bees(data), p in plants(data)]) ) / prod(length.([bees(data), plants(data)]))
=#
