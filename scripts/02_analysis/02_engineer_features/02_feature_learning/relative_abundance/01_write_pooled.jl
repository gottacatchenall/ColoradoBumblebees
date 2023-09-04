using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()
    emb = representations(data, Pooled())

    ColoradoBumblebees.save(emb)
end 

main()