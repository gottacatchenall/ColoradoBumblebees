using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    for s in species 
        sdms = make_sdms(s, occurrence_df; cluster=cluster)
    
        for (i,sdm) in enumerate(sdms)
            sdm_dir = sdmdir(sdm)
            mkpath(sdm_dir)
            @info sdm_dir
            SpeciesDistributionToolkit.save(joinpath(sdm_dir, "prediction.tif"), sdm.probability)
            SpeciesDistributionToolkit.save(joinpath(sdm_dir, "uncertainty.tif"), sdm.uncertainty)
            json_string = JSON.json(sdm.fit_stats)
            open(joinpath(sdm_dir, "fit.json"), "w") do f
                JSON.print(f, json_string)
            end
            #SpeciesDistributionToolkit.save("$(s)_tmp$i.tif", sdm.probability)
            #ColoradoBumblebees.save(sdm)
        end
    end 

end

main()
