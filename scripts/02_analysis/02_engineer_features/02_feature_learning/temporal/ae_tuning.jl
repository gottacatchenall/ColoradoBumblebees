using DrWatson
@quickactivate :ColoradoBumblebees

using CUDA

function main(;unit=RNN)
    print("GPU: $(ColoradoBumblebees.GPU_AVAILABLE)")
    data = load_data()

    d = Dict(
        :rnn_dims => [ [1, 8, 1], [1, 16, 1], [1, 32, 1], 
                       [1, 8, 8, 1], [1, 16, 16, 1], [1, 32, 32, 1], 
                       [1, 8, 8, 8, 1], [1, 16, 16, 16, 1], [1, 32, 32, 32, 1]],
        :encoder_dims => [[TEMPORAL_INPUT_DIM, 4], [TEMPORAL_INPUT_DIM, 8], [TEMPORAL_INPUT_DIM, 16]],
        :decoder_dims => [[4, TEMPORAL_INPUT_DIM], [8, TEMPORAL_INPUT_DIM], [16, TEMPORAL_INPUT_DIM]],
    )

    a = dict_list(d)
    
    enc_dec_sets = filter(x-> x[:encoder_dims][end] == x[:decoder_dims][begin], a)


    models = vcat([RecurrentAutoencoder{Standard}(;  unit=unit,opt=ADAM(0.005), e...) for e in enc_dec_sets],
    [RecurrentAutoencoder{Variational}(; unit=unit, opt=ADAM(0.001), e...) for e in enc_dec_sets])


    for m in models
        rep = representations(data, m)
        ColoradoBumblebees.save(rep)
    end
end 


main(unit=LSTM)
#main(unit=RNN)
#main(unit=GRU)