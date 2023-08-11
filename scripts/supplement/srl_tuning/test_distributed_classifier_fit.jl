# Defining variables available (and independent) across all processes

@everywhere begin
    using DrWatson
end

@everywhere begin
    @quickactivate :ColoradoBumblebees
    data = load_data()
end 

@everywhere begin  
    
    function run_model(embed; num_replicates = 64)
        model = RandomForest()
        feat_df = feature_dataframe(data, embed)
        bf = batch_fit(model, embed, feat_df, num_replicates)
        @info "Fit done"
        ColoradoBumblebees.save(bf)
    end

end

# Main script
using SlurmClusterManager
addprocs(SlurmManager())

@everywhere println("hello from $(myid()):$(gethostname())")


#=
srdir = joinpath(artifactdir(), "species_representations")
reps = readdir(srdir)


status = @showprogress pmap(1:2) do i
    try
        run_model(ColoradoBumblebees.load(joinpath(srdir, reps[i])))
        true # success
    catch e
        @info e
        false # failure
    end
end=#
  
  