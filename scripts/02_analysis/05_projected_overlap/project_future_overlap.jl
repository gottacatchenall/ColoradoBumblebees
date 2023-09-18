using DrWatson
@quickactivate :ColoradoBumblebees

function main()    
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

    SSPs = [SSP1_26, SSP2_45, SSP3_70]
    all_treatments = vcat(Dict(:ssp=>Baseline, :time=>baseline()), dict_list(Dict(:ssp=>SSPs, :time=>TIMESPANS[2:end]))) 
    this_ssp, this_timespan  = all_treatments[job_id][:ssp], all_treatments[job_id][:time]

    treatment_str = string(this_ssp)*"_"*string(this_timespan)

    # Saving using ArchGDAL sometimes writes scratch files to local dir during
    # the process of writing. Idk if its due to lack of RAM or something else,
    # but running each species in its own dir stops simulataneous overwriting of
    # scratch files
    mkpath("./tmp/$(treatment_str)")
    cd("./tmp/$(treatment_str)")
    @info "Job: $job_id, Species: $(treatment_str), WD: $(pwd())"

    overlap = compute_overlap(BEST_FIT_DIR, this_timespan, this_ssp)

    ColoradoBumblebees.save(overlap; cluster=ColoradoBumblebees.CLUSTER)
end
 
main()