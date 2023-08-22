using DrWatson
@quickactivate :ColoradoBumblebees

function main()
    data = load_data()
    km = KMeansEnvironmentEmbedding
    mods = [km(i) for i in 1:10]
    embs = [representations(data, m) for m in mods]
    for e in embs
        ColoradoBumblebees.save(e)
    end
end

main()