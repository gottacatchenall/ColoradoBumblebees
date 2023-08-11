using DrWatson
@quickactivate :ColoradoBumblebees


function run_model()
    rep_dir = joinpath(artifactdir(), "species_representations")

    temporal_reps = sort(filter(x->contains(x, "RecurrentAutoencoder"), readdir(rep_dir)))
    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

    jp = joinpath(rep_dir, temporal_reps[job_id])
    this_embed = ColoradoBumblebees.load(jp)

    data = load_data()
    model = RandomForest()
    feat_df = feature_dataframe(data, this_embed)
    bf = batch_fit(model, this_embed, feat_df, num_replicates)
    ColoradoBumblebees.save(bf)
end


run_model()