using DrWatson
@quickactivate :ColoradoBumblebees

function load_baseline_overlap(cluster)
    lead = cluster ? "/scratch/mcatchen/artifacts/projected_overlap" : joinpath(artifactdir(), "projected_overlap")
    baseline_path = filter(x->contains(x, "Baseline"), readdir(lead))[1]
    ColoradoBumblebees.load(joinpath(lead, baseline_path))
end

function main()
    best_rep = "RelativeAbundance_Structural_Environment"
    best_fit_dir = joinpath(artifactdir(), "classification_fits", "multiple_representations", "brt", best_rep)

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

    SSPs = [SSP1_26, SSP2_45, SSP3_70]
    all_treatments = dict_list(Dict(:ssp=>SSPs, :time=>TIMESPANS[2:end]))
    
    this_ssp, this_timespan  = all_treatments[job_id][:ssp], all_treatments[job_id][:time]

    treatment_str = string(this_ssp)*"_"*string(this_timespan)

    # Saving using ArchGDAL sometimes writes scratch files to local dir during
    # the process of writing. Idk if its due to lack of RAM or something else,
    # but running each species in its own dir stops simulataneous overwriting of
    # scratch files
    mkpath("./tmp/$(treatment_str)")
    cd("./tmp/$(treatment_str)")
    @info "Job: $job_id, Species: $(treatment_str), WD: $(pwd())"

    baseline_po = load_baseline_overlap(ColoradoBumblebees.CLUSTER)
    projected_overlap = compute_future_overlap(baseline_po, best_fit_dir, this_timespan, this_ssp)

    ColoradoBumblebees.save(projected_overlap; cluster=ColoradoBumblebees.CLUSTER)
end
 
main()