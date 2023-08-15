using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()

    param_dict = Dict(
        :embedding_dim => [4:4:16],
        :number_of_walks => [50, 100, 250, 500],
        :walk_length => [10, 50, 100]
    )

    treatments = dict_list(param_dict)


    embs = [representations(data, PhylogeneticNode2Vec())]

    for e in embs
        ColoradoBumblebees.save(e)
    end
end 

data = load_data()
x = representations(data, PhylogeneticNode2Vec(embedding_dim=16))
bf = batch_fit(RandomForest(), x, feature_dataframe(data, x), 64)


praucs(bf)