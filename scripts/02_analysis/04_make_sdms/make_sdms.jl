using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    for s in species 
        sdms = make_sdms(s, occurrence_df; cluster=cluster)
    
        for (i,sdm) in enumerate(sdms)
            SpeciesDistributionToolkit.save("$(s)_tmp$i.tif", sdm.probability)
            #ColoradoBumblebees.save(sdm)
        end
    end 

end

main()
