using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()

    param_dict = Dict(
        :numtraits => [100, 1000, 10000], 
        :stddev_of_stddev => [1.0,2.0,3.0],
        :truncated_dims => [i for i in 2:2:16]
    )

    treatments = dict_list(param_dict)


    embs = [representations(data, SimulatedTraits(numtraits=r[:numtraits], stddev_of_stddev=r[:stddev_of_stddev], truncated_dims=r[:truncated_dims])) for r in treatments]

    for e in embs
        ColoradoBumblebees.save(e)
    end
end 

main()