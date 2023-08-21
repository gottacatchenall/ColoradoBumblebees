using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()

    param_dict = Dict(
        :numtraits => [100, 1000, 10000], 
        :variance_distribution => [TruncatedNormal(0.,σ, 0, Inf) for σ in [1.0,2.0,3.0]],
        :truncated_dims => [i for i in 2:2:16]
    )

    treatments = dict_list(param_dict)


    embs = [representations(data, SimulatedTraits(numtraits=r[:numtraits], variance_distribution=r[:variance_distribution], truncated_dims=r[:truncated_dims])) for r in treatments]

    for e in embs
        ColoradoBumblebees.save(e)
    end
end 

data = load_data()
x = representations(data, SimulatedTraits(500, TruncatedNormal(0,1,0,Inf), 16))

bf = batch_fit(RandomForest(), x, feature_dataframe(data, x), 64)

praucs(bf)
