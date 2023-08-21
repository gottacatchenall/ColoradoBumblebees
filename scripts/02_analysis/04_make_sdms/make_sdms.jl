using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_species = species[job_id]

    mkpath("./tmp/$(this_species)")
    cd("./tmp/$(this_species)")
    @info "Job: $job_id, Species: $(this_species), WD: $(pwd())"


    sdms = make_sdms(this_species, occurrence_df; cluster=cluster, pa_buffer_distance=8)
    for (i,sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        

        # For some reason theres weird saving issues. 
        # Attempt 1: assume it is a GDAL <-> filesystem fuckup, and try writing
        # the same sdm a few times.

        # They all write prediction local scratch files, maybe each job should
        # cd into its own species specific dir?
        
        succeeded = false
        for attempt in 1:50
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
