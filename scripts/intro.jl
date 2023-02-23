using DrWatson
@quickactivate "ColoradoBumblebees"

using MLJ
using ProgressMeter
using DataFrames

# Here you may include files from the source directory
include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()


embeddings = [KMeansEnvironmentEmbedding(15)]

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
ens_size = 100

Is = shuffle(1:nrow(df))
cut = Int32(floor(0.7*nrow(df)))
Itrain, Itest  = Is[1:cut], Is[cut+1:end]


ytest = [x == true for x in y[Itest]]

rf = RandomForest()
dt = DecisionTree()
brt = BRT()

ypredict = zeros(length(Itest))
models = []
@showprogress for i in 1:ens_size
    mach = machine(dt,X,y)
    theserows = balance_sample(y, Itrain, 128, 0.5)
    fit!(mach, rows=theserows, verbosity=0)
    pred = predict(mach, rows=Itest)

    ypredict .+= [p.prob_given_ref[2] for p in pred]
    #push!(models, mach)
end

ypredict = ypredict ./ ens_size

computemeasures(ytest, ypredict)