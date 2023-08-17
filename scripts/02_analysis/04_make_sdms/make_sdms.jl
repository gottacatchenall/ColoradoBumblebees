using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    for s in species 
        sdms = make_sdms(s, occurrence_df; cluster=cluster)
        
        SpeciesDistributionToolkit.save("tmp.tif", sdms[1].probability)
        return 

        for sdm in sdms
            ColoradoBumblebees.save(sdm)
        end
    end 

end

main()
