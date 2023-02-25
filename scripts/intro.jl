using DrWatson
@quickactivate "ColoradoBumblebees"

using MLJ
using ProgressMeter
using DataFrames
using Flux
using Distributions

# Here you may include files from the source directory
include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()


embeddings = [
    KMeansSpatialEmbedding(5), 
    Autoencoder{Variational}(
        unit=LSTM, 
        η=5e-3, 
        n_epochs=250, 
        encoder_dims=[TEMPORAL_INPUT_DIM,256,64,16,4], 
        decoder_dims=[4,64,TEMPORAL_INPUT_DIM], 
        train_proportion=1.),

    Autoencoder{Standard}(
        unit=Dense, 
        η=5e-3, 
        n_epochs=2000,
        encoder_dims=[TEMPORAL_INPUT_DIM,256,64,16,4], 
        decoder_dims=[4,64,TEMPORAL_INPUT_DIM],
        train_proportion=1.),
    KMeansEnvironmentEmbedding(5),
]

df = feature_dataframe(data, embeddings)

y, X, _, = unpack(df, ==(:interaction), ∉([:bee,:plant]); rng = 123)
y = coerce(y, Multiclass{2})

DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
BRT = @load EvoTreeClassifier pkg = EvoTrees

rf = RandomForest()
dt = DecisionTree()
brt = BRT()


balances = 0.2:0.1:8
batch_sizes = [2^i for i in 3:9]
ensemble_sizes = [2^i for i in 3:9]
reps = 64

results = [DataFrame(prauc=[], rocauc=[], balance=[], batch_size=[], ensemble_size=[]) for i in 1:reps]

Threads.@threads for r in 1:reps
    Is = shuffle(1:nrow(df))
    cut = Int32(floor(0.8*nrow(df)))
    Itrain, Itest  = Is[1:cut], Is[cut+1:end]
    ytest = [x == true for x in y[Itest]]

    for b in balances
        for bsize in batch_sizes
            for es in ensemble_sizes          
                @info "Balance: $(b)\tBatch Size: $(bsize)\tEnsemble Size$(es)\tThread:$(Threads.threadid())"     
                ypredict = zeros(length(Itest))
                for i in 1:es
                    mach = machine(dt,X,y)
                    theserows = balance_sample(y, Itrain, bsize, b)
                    fit!(mach, rows=theserows, verbosity=0)
                    pred = predict(mach, rows=Itest)
                    ypredict .+= [p.prob_given_ref[2] for p in pred]            
                end

                ypredict = ypredict ./ es
                meas = computemeasures(ytest, ypredict)
                push!(results[Threads.threadid()].prauc, meas[:prauc])
                push!(results[Threads.threadid()].rocauc, meas[:rocauc])
                push!(results[Threads.threadid()].balance, b)
                push!(results[Threads.threadid()].ensemble_size, es)
                push!(results[Threads.threadid()].batch_size, bsize)
            end
        end
    end
end


# Simulation suffix
_jobid = get(ENV, "SLURM_ARRAY_TASK_ID", 1)
_jobcount = get(ENV, "SLURM_ARRAY_TASK_COUNT", 1)

CSV.write(datadir("output_$(_jobid).csv"), vcat(results...))


#=
Is = shuffle(1:nrow(df))
cut = Int32(floor(0.8*nrow(df)))
Itrain, Itest  = Is[1:cut], Is[cut+1:end]
ytest = [x == true for x in y[Itest]]
ens_size = 256
ypredict = zeros(length(Itest))
for i in 1:ens_size
    for m in models
        mach = machine(m,X,y)
        theserows = balance_sample(y, Itrain, 64, 0.5)
        fit!(mach, rows=theserows, verbosity=0)
        pred = predict(mach, rows=Itest)
        ypredict .+= [p.prob_given_ref[2] for p in pred]
    end 
end


ypredict = ypredict ./ (length(models)*ens_size)
computemeasures(ytest, ypredict) 

=#