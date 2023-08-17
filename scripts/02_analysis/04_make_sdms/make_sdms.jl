using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    for s in species 
        sdms = make_sdms(s, occurrence_df; cluster=cluster)
        
        return 

        for (i,sdm) in enumerate(sdms)
            SpeciesDistributionToolkit.save("$s_tmp$i.tif", sdm.probability)
            #ColoradoBumblebees.save(sdm)
        end
    end 

end

main()
