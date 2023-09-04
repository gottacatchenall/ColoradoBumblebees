using DrWatson
@quickactivate :ColoradoBumblebees


function main(num_replicates)
    treatments = filter(x->length(x)>0,collect(powerset(BEST_REPRESENTATIONS)))
    treatments = treatments[sortperm([string(t) for t in treatments])]
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_treatment = treatments[job_id]

    model = RandomForest()
    bf = compare_representations(model, num_replicates, this_treatment)
    @info "Saving to $(path(bf))"
    ColoradoBumblebees.save(bf)
end 

main(128) # 1371 seconds