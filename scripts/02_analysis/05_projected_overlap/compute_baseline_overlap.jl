using DrWatson
@quickactivate :ColoradoBumblebees

function main()
    best_rep = "RelativeAbundance_Structural_Environment"
    best_fit_dir = joinpath(artifactdir(), "classification_fits", "multiple_representations", "brt", best_rep)
    baseline_po = compute_overlap(best_fit_dir, baseline(), Baseline)

    ColoradoBumblebees.save(baseline_po; cluster=ColoradoBumblebees.CLUSTER)
end 

