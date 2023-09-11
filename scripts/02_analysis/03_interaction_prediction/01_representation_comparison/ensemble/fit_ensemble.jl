using DrWatson
@quickactivate :ColoradoBumblebees

function main(num_replicates; outdir="")
    treatments = filter(x->length(x)>0,collect(powerset(BEST_REPRESENTATIONS)))
    treatments = treatments[sortperm([string(t) for t in treatments])]
    
    n = [join(string.(supertype.(typeof.(t))), "_") for t in treatments]    
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_treatment = treatments[job_id]
    this_name = n[job_id]
    this_outpath = joinpath(outdir, "classification_fits", "multiple_representations", "ensemble", this_name)


    bf = compare_representation_ensemble(num_replicates, this_treatment)
    ColoradoBumblebees.save(bf; outdir=this_outpath)
end 

main(128; outdir="/scratch/mcatchen") 