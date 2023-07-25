using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

using Distributions
using Flux
using Dates
using DataFrames
using Random
using MLJ
using ProgressMeter
using CSV
using CairoMakie

const RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
const XGBoostClassifier = @load XGBoostClassifier verbosity = 0



data = load_data()

ra = RecurrentAutoencoder(
    rnn_dims=[1,8,1], 
    encoder_dims=[TEMPORAL_INPUT_DIM,8], 
    decoder_dims=[8,TEMPORAL_INPUT_DIM], 
    unit=LSTM)

df = feature_dataframe(data, [ra])


df = feature_dataframe(data, [LFSVD()])


y, X, species_pairs = unpack(df, ==(:interaction), ∉([:bee, :plant]); rng=123)

X = MLJ.transform(fit!(machine(Standardizer(), X)),X)


y = coerce(y, Multiclass{2})

rf = RandomForest()

mach, cv = crossvalidation(rf, X, y, species_pairs)