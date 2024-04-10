using DrWatson
@quickactivate :ColoradoBumblebees

const RADIUS = 400.0
const BUFFER = 20.0

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


    sdms = make_sdms(
        this_species,
        occurrence_df; 
        cluster=cluster, 
        thickening_distance=RADIUS, 
        buffer=BUFFER
    )
    
    # sometimes the default buffer, thickening radius, and bias give weird fits.
    # for the sake of brevity, lets hope its a few species with similar
    # properties
    if sdms[1].fit_stats[:rocauc] < 0.5  
        iostream =  Base.open("../../sdms_to_fix.txt","a")
        write(iostream, "$(s.name)\n");
        return 
    end 

    for (i, sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        ColoradoBumblebees.save(sdm)
    end
end

main()


