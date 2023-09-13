using DrWatson
@quickactivate :ColoradoBumblebees


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


    sdms = make_sdms(this_species, occurrence_df; cluster=cluster, pa_buffer_distance=8)
    for (i,sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        ColoradoBumblebees.save(sdm)
    end


end
 
#main()



#=
 [ Info: Bombus californicus is broken
[ Info: Dodecatheon pulchellum is broken
[ Info: Ericameria parryi is broken
[ Info: Iris missouriensis is broken
[ Info: Pedicularis bracteosa is broken
[ Info: Verbascum thapsus is broken
=#

function main2()

    spnames = ["Bombus californicus", "Dodecatheon pulchellum", "Ericameria parryi", "Iris missouriensis", "Pedicularis bracteosa", "Verbascum thapsus"]

    for sp in spnames
        sdms = make_sdms(sp, load_occurrence_data(); cluster=false, pa_buffer_distance=8)
        for (i,sdm) in enumerate(sdms)
            sdm_dir = sdmdir(sdm)
            mkpath(sdm_dir)
            ColoradoBumblebees.save(sdm)
        end
    end 

end 

main2()


