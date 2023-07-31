using DrWatson
@quickactivate :ColoradoBumblebees

const RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
const XGBoostClassifier = @load XGBoostClassifier verbosity = 0

data = load_data()


d = Dict(
    :rnn_dims => [[1, 8, 1], [1, 16, 1], [1, 32, 1], [1, 8, 8, 1], [1, 16, 16, 1], [1, 32, 32, 1]],
    :encoder_dims => [[TEMPORAL_INPUT_DIM, 4], [TEMPORAL_INPUT_DIM, 8], [TEMPORAL_INPUT_DIM, 16]],
    :decoder_dims => [[4, TEMPORAL_INPUT_DIM], [8, TEMPORAL_INPUT_DIM], [16, TEMPORAL_INPUT_DIM]],
)

a = dict_list(d)

filter(x-> x[:encoder_dims][end] == x[:decoder_dims][begin], a)





# Now filter out the sets where the encoder and decoder dimensions don't match.

ra = RecurrentAutoencoder(
    rnn_dims=[1,8,1], 
    encoder_dims=[TEMPORAL_INPUT_DIM,8], 
    decoder_dims=[8,TEMPORAL_INPUT_DIM], 
    unit=LSTM)

df = feature_dataframe(data, [ra])




df = feature_dataframe(data, [MetawebSVD(truncation_dims=4, embed_dims=4)])


y, X, species_pairs = unpack(df, ==(:interaction), ∉([:bee, :plant]); rng=123)

X = MLJ.transform(fit!(machine(Standardizer(), X)),X)


y = coerce(y, Multiclass{2})

rf = RandomForest()


mach, cv = crossvalidation(rf, X, y, species_pairs)