function _get_srl_dir(emb)
    "SpeciesRepresenetation_$(typeof(emb))_"*savename(emb)
end 

# Note: This function expects you to be ussing SLURM array lobs
function compare_representations(model, num_replicates, this_treatment)

    rep_dir = joinpath(artifactdir(), "species_representations")

    reps = [ColoradoBumblebees.load(joinpath(rep_dir, _get_srl_dir(emb))) for emb in this_treatment]

    data = load_data()
    feat_df = feature_dataframe(data, reps)

    bf = batch_fit(model, reps, feat_df, num_replicates)
    return bf
end


# Best representations. This will probably become a constant?
# Phylogenetic: 
    # Node2Vec: Walks: 250, Walklength: 100, 16: PRAUC: 0.63
    # Simulated traits: Trait σ: 1, Embed Dims: 14, PRAUC: 0.61
# Structural 
    # SVD: Trunc 14, Embed 2: PRAUC 0.818
    # LFSVD: 0.815, Trunc 6, Embed 6
# Environment 
    # K = 7, PR 0.63
# Temporal 
    # Standard GRU 
    # 1 hidden layer, 8 channels, 16 embed,  PR 0.77
# Relative Abd: Pooled
