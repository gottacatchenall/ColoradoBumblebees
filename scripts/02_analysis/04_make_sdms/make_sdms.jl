using DrWatson
@quickactivate :ColoradoBumblebees

const RADIUS = 400.0
const BUFFER = 40.0

function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_species = species[job_id]

    # Saving using ArchGDAL sometimes writes scratch files to local dir during
    # the process of writing. Idk if its due to lack of RAM or something else,
    # but running each species in its own dir stops simulataneous overwriting of
    # scratch files
    mkpath("./tmp/$(this_species)")
    cd("./tmp/$(this_species)")
    @info "Job: $job_id, Species: $(this_species), WD: $(pwd())"


    sdms = make_sdms(this_species, occurrence_df; cluster=cluster, radius=RADIUS, buffer=BUFFER)
    for (i, sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        ColoradoBumblebees.save(sdm)
    end
end

main()


