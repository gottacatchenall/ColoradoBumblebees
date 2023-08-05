using DrWatson
@quickactivate :ColoradoBumblebees


data = load_data()

d = Dict(
    :rnn_dims => [[1, 8, 1], [1, 16, 1], [1, 32, 1], [1, 8, 8, 1], [1, 16, 16, 1], [1, 32, 32, 1]],
    :encoder_dims => [[TEMPORAL_INPUT_DIM, 4], [TEMPORAL_INPUT_DIM, 8], [TEMPORAL_INPUT_DIM, 16]],
    :decoder_dims => [[4, TEMPORAL_INPUT_DIM], [8, TEMPORAL_INPUT_DIM], [16, TEMPORAL_INPUT_DIM]],
)

a = dict_list(d)
filter(x-> x[:encoder_dims][end] == x[:decoder_dims][begin], a)


rep.embedding

lfsvd = representations(data, LFSVD(embed_dims=4))
phy = representations(data, SimulatedTraits())
env = representations(data, KMeansEnvironmentEmbedding())

feat_df = feature_dataframe(data, lfsvd)
xgb = XGBoost()
@time fit_classifier(xgb, lfsvd, feat_df; train_proportion=0.7)

# Test save and load of representations
ColoradoBumblebees.save(lfsvd)
dir = path(lfsvd)
ColoradoBumblebees.load(dir)



# Test save and load of BatchFit

# Single Rep
reps =  [lfsvd, phy, env]
feat_df = feature_dataframe(data, lfsvd)
bf = @time batch_fit(xgb, lfsvd, feat_df, 8)

save(bf)
dir = path(bf)
ColoradoBumblebees.load(dir)


# Many reps
feat_df = feature_dataframe(data, reps)
bf = @time batch_fit(xgb, reps, feat_df, 8)

save(bf)
dir = path(bf)
ColoradoBumblebees.load(dir)












# Now filter out the sets where the encoder and decoder dimensions don't match.

 ra = RecurrentAutoencoder(
    rnn_dims=[1,8,1], 
    encoder_dims=[TEMPORAL_INPUT_DIM,8], 
    decoder_dims=[8,TEMPORAL_INPUT_DIM], 
    unit=LSTM)

df = feature_dataframe(data, [ra])


representations(data, MetawebSVD())


df = feature_dataframe(data, [MetawebSVD(truncation_dims=4, embed_dims=4)])


y, X, species_pairs = unpack(df, ==(:interaction), ∉([:bee, :plant]); rng=123)

X = MLJ.transform(fit!(machine(Standardizer(), X)),X)


y = coerce(y, Multiclass{2})

mach = machine(RandomForest()(), X,y)

fit!(mach)

y_predict = MLJ.predict(mach)


rf = RandomForest()

prediction, cv = fit_model(df, rf)