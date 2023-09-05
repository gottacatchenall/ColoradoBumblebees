using DrWatson
@quickactivate :ColoradoBumblebees

function main(num_replicates)
    treatments = filter(x->length(x)>0,collect(powerset(BEST_REPRESENTATIONS)))
    treatments = treatments[sortperm([string(t) for t in treatments])]
    names = [join(string.(supertype.(typeof.(t))), "_") for t in treatments]    
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_treatment = treatments[job_id]
    this_name = names[job_id]
    this_outpath = joinpath(artifactdir(), "classification_fits", "multiple_representations", "logistic", this_name)


    model = Logistic()
    bf = compare_representations(model, num_replicates, this_treatment)
    ColoradoBumblebees.save(bf; outpath = this_outpath)
end 

main(128) 