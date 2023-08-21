using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees
using MLJ
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

function _make_feature_df(feat)
    labeldf = label_dataframe(data)

    feat_dim = length(values(feat) |> first)
    cols = []
    push!(
        cols,
        [
            "BEE$fi" => zeros(nrow(labeldf)) for fi in 1:feat_dim
        ]...,
    )
    push!(
        cols,
        [
            "PLANT$fi" => zeros(nrow(labeldf)) for fi in 1:feat_dim
        ]...,
    )

    featdf = DataFrame(cols...)

    for (ri, r) in enumerate(eachrow(labeldf))
        b, p = r.bee, r.plant
        beevec, plantvec = feat[b], feat[p]
        for i in 1:feat_dim
            featdf[ri, "BEE$i"] = beevec[i]
        end
        for i in 1:feat_dim
            featdf[ri, "PLANT$i"] = plantvec[i]
        end
    end
    return hcat(labeldf, featdf)
end 

function get_feature_dict(result_df)
    dict = Dict()
    data = load_data()
    for r in eachrow(result_df)
        sp = r.species
        spobj = contains(sp, "Bombus") ? bee(data, sp) : plant(data, sp)
        merge!(dict, Dict(spobj=>Array(r[1:end-1])))
    end
    return dict 
end 

function get_feature_df(result_df)
    _make_feature_df(get_feature_dict(result_df))
end 

function fit(feature_df;n_reps=16)
    y, X, species_pairs = unpack(feature_df, ==(:interaction), ∉([:bee, :plant]); rng=123)

    X = MLJ.transform(fit!(machine(Standardizer(), X)),X)
    y = coerce(y, Multiclass{2})

    s = 0
    for i in 1:n_reps
        fitdf, _ = crossvalidation(rf, X, y, species_pairs)
        s +=  fitdf.prauc[1]
    end 
    @info "Mean PRAUC: $(s/n_reps)"
end


result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,8,1)_ENC(147,16)_DEC(16,147).csv" ,DataFrame) 
# Mean PRAUC: 0,7414 
fit(get_feature_df(result_df); n_reps=64)

result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,32,1)_ENC(147,16)_DEC(16,147).csv" ,DataFrame) 
# Mean PRAUC: 0.7306
fit(get_feature_df(result_df))

result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,8,1)_ENC(147,32)_DEC(32,147).csv" ,DataFrame) 
# Mean PRAUC: 0.7476
fit(get_feature_df(result_df))


result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,8,8,1)_ENC(147,16)_DEC(16,147).csv" ,DataFrame) 
# Mean PRAUC: 0.745
fit(get_feature_df(result_df), n_reps=64)


result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,8,8,8,1)_ENC(147,16)_DEC(16,147).csv" ,DataFrame) 
# Mean PRAUC: 0.732
fit(get_feature_df(result_df), n_reps=64)


result_df = CSV.read("scripts/supplement/temporal_autoencoder_model_tuning/LSTM(1,32,32,1)_ENC(147,16)_DEC(16,147).csv" ,DataFrame) 
# Mean PRAUC: 0.739
fit(get_feature_df(result_df); n_reps=64)