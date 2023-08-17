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
        ColoradoBumblebees.save(sdm)
    end

end

main()
