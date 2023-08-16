using DrWatson
@quickactivate :ColoradoBumblebees


function main()
    data = load_data()
    occurrence_df = load_occurrence_data()

    species = sort(unique(occurrence_df.species))
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

    @info "Species: $(species[job_id])"
    cluster = ColoradoBumblebees.CLUSTER

    lk = ReentrantLock()

    sdms = make_sdms(species[job_id], occurrence_df; cluster=cluster)
    
    lock(lk) do
        for sdm in sdms
            ColoradoBumblebees.save(sdm)
        end
    end 

end

main()
