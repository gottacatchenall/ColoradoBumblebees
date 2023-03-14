using DrWatson
@quickactivate "ColoradoBumblebees"

using MLJ
using ProgressMeter
using DataFrames
using Flux
using Distributions
using ParameterSchedulers
using ParameterSchedulers: Scheduler
using CSV

# Here you may include files from the source directory
include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

embeddings = [
    SimulatedTraits(),
 #   Pooled(),
    # Hierarchical(),
    MetawebSVD(8)
 #   KMeansSpatialEmbedding(3),
    # KMeansEnvironmentEmbedding(3),

    #= Autoencoder{Variational}(
         unit=LSTM, 
         n_epochs=250, 
         dropout=0.,
         opt=ADAM(0.001),
         encoder_dims=[TEMPORAL_INPUT_DIM,128,64,8], 
         decoder_dims=[8,64,128,TEMPORAL_INPUT_DIM], 
         train_proportion=1.)

     Autoencoder{Standard}(
         unit=Dense, 
         n_epochs=2500,
         opt=ADAM(1e-4),
         encoder_dims=[TEMPORAL_INPUT_DIM,32,4], 
         decoder_dims=[4,32,TEMPORAL_INPUT_DIM],
         train_proportion=1.) =#
]

df = feature_dataframe(data, embeddings)

y, X, _, = unpack(df, ==(:interaction), ∉([:bee, :plant]); rng=123)
y = coerce(y, Multiclass{2})

single_run(X, y, 256, 128)

DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
BRT = @load EvoTreeClassifier pkg = EvoTrees

# brt = BRT(n_subfeatures=20, n_trees=100, max_depth=30, nbins=64)
brt = BRT(; n_subfeatures=20)
rf = RandomForest(; n_subfeatures=20, n_trees=100)
#
MLJ.evaluate(rf, X, y; measure=computemeasures_mlj)

function single_run(X, y, ens_size=256, batch_size=64)
    rf = RandomForest()
    Is = shuffle(1:nrow(df))
    cut = Int32(floor(0.8 * nrow(df)))
    Itrain, Itest = Is[1:cut], Is[(cut + 1):end]
    ytest = [x == true for x in y[Itest]]
    ypredict = zeros(length(Itest))
    for i in 1:ens_size
        mach = machine(rf, X, y)
        theserows = balance_sample(y, Itrain, batch_size, 0.5)
        fit!(mach; rows=theserows, verbosity=0)
        pred = MLJ.predict(mach; rows=Itest)
        ypredict .+= [p.prob_given_ref[2] for p in pred]
    end

    ypredict = ypredict ./ (ens_size)
    return computemeasures(ytest, ypredict)
end

balances = 0.2:0.05:0.8
batch_sizes = [2^i for i in 3:7]
ensemble_sizes = [2^i for i in 3:9]
reps = Threads.nthreads()

results = [
    DataFrame(; prauc=[], rocauc=[], balance=[], batch_size=[], ensemble_size=[]) for
    i in 1:reps
]

Threads.@threads for r in 1:reps
    Is = shuffle(1:nrow(df))
    cut = Int32(floor(0.8 * nrow(df)))
    Itrain, Itest = Is[1:cut], Is[(cut + 1):end]
    ytest = [x == true for x in y[Itest]]

    for b in balances
        for bsize in batch_sizes
            for es in ensemble_sizes
                @info "Balance: $(b)\tBatch Size: $(bsize)\tEnsemble Size$(es)\tThread:$(Threads.threadid())"
                ypredict = zeros(length(Itest))
                for i in 1:es
                    mach = machine(dt, X, y)
                    theserows = balance_sample(y, Itrain, bsize, b)
                    fit!(mach; rows=theserows, verbosity=0)
                    pred = predict(mach; rows=Itest)
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
