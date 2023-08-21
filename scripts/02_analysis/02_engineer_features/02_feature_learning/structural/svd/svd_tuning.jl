using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()

    param_dict = Dict(
        :model => [MetawebSVD, LFSVD], 
        :truncate_dims => [i for i in 2:4:18],
        :embed_dims => [i for i in 2:4:18]
    )
    treatments = filter(x-> x[:truncate_dims] >= x[:embed_dims], dict_list(param_dict))
    models =  [r[:model](embed_dims=r[:embed_dims], truncation_dims=r[:truncate_dims]) for r in treatments]

    embs = [representations(data, m) for m in models]

    for e in embs
        ColoradoBumblebees.save(e)
    end 
end

main()