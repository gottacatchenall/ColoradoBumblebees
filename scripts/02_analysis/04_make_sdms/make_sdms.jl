using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_species = species[job_id]
    @info "Job: $job_id, Species: $(this_species)"

    sdms = make_sdms(this_species, occurrence_df; cluster=cluster)
    for (i,sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        
        # For some reason theres weird saving issues. 
        # Attempt 1: assume it is a GDAL <-> filesystem fuckup, and try writing
        # the same sdm a few times.
        
        # This is not ideal, I still think it might be the `similar` call when
        # constructing the bg point mask. idk

        succeeded = false
        for attempt in 1:10
            if succeeded
                break
            end
            try 
                ColoradoBumblebees.save(sdm)
                succeeded = true
            catch e
                @info "Failed on $(this_species) attempt $attempt"
            end
        end 
    end


end

main()
