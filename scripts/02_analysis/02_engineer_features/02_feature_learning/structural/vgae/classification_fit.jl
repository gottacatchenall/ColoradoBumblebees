using DrWatson
@quickactivate :ColoradoBumblebees

function main()
    rep_dir = joinpath(artifactdir(), "species_representations")
    gae_reps = sort(filter(x->contains(x, "GraphAutoencoder"), readdir(rep_dir)))
   
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

    jp = joinpath(rep_dir, gae_reps[job_id])
    this_rep = ColoradoBumblebees.load(jp)

    data = load_data()
    feat_df = feature_dataframe(data, this_rep)

    model = RandomForest()

    bf = batch_fit(model, this_rep, feat_df, num_replicates)
    ColoradoBumblebees.save(bf)
end