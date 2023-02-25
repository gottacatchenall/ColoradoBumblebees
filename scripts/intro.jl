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
    Autoencoder{Variational}(unit=LSTM, η=1e-3, n_epochs=250, encoder_dims=[TEMPORAL_INPUT_DIM,256,64,16,8], decoder_dims=[8,64,TEMPORAL_INPUT_DIM], train_proportion=1.),
    Autoencoder{Standard}(unit=Dense, η=1e-2, encoder_dims=[TEMPORAL_INPUT_DIM,256,64,16,8], decoder_dims=[8,64,TEMPORAL_INPUT_DIM], train_proportion=1.),
    #    KMeansEnvironmentEmbedding(10)
    #       PCA them instead of whitening
]

#embeddings = [KMeansSpatialEmbedding(20)]

df = feature_dataframe(data, embeddings)


y, X, _, = unpack(df, ==(:interaction), ∉([:bee,:plant]); rng = 123)
y = coerce(y, Multiclass{2})

DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
BRT = @load EvoTreeClassifier pkg = EvoTrees

function balance_sample(y, I, batch_size=64, true_pct=0.5)
    itrue = findall(x-> x == true, y)
    ifalse = findall(x-> x == false, y)

    filter!(x->x ∈ I, itrue)
    filter!(x->x ∈ I, ifalse)

    ntrue, nfalse = Int32(floor(batch_size*true_pct)),Int32(floor(batch_size*(1-true_pct)))
    [shuffle(itrue)[1:ntrue]..., shuffle(ifalse)[1:nfalse]...]
end

## testing method for balanced training 

ytest = [x == true for x in y[Itest]]


SVM = @load SVMClassifier pkg=ScikitLearn

rf = RandomForest()
dt = DecisionTree()
brt = BRT()
svm = SVM()

models = [dt]


res_df = DataFrame(prauc=[], rocauc=[], balance=[], batch_size=[], ensemble_size=[])

balances = 0.2:0.1:8
batch_sizes = [16,64,128]
ensemble_sizes = [64,128,256]
reps = 10

@showprogress for b in balances
    for bsize in batch_sizes
        for es in ensemble_sizes
            for r in 1:reps
                @info (b,bsize,es,r)
                Is = shuffle(1:nrow(df))
                cut = Int32(floor(0.8*nrow(df)))
                Itrain, Itest  = Is[1:cut], Is[cut+1:end]

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
                push!(res_df.prauc, meas[:prauc])
                push!(res_df.rocauc, meas[:rocauc])
                push!(res_df.balance, b)
                push!(res_df.ensemble_size, es)
                push!(res_df.batch_size, bsize)

            end
        end
    end
end






ens_size = 128
ypredict = zeros(length(Itest))
 for i in 1:ens_size
    for m in models
        mach = machine(m,X,y)
        theserows = balance_sample(y, Itrain, 16, 0.25)
        fit!(mach, rows=theserows, verbosity=0)
        pred = predict(mach, rows=Itest)
        ypredict .+= [p.prob_given_ref[2] for p in pred]
    end 
end


ypredict = ypredict ./ (length(models)*ens_size)
computemeasures(ytest, ypredict) 

