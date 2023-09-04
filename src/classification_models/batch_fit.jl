#=
function get_metadata(mach::Machine)
    fn = fieldnames(typeof(mach.model))
    dict = Dict()
    for f in fn
        dict = merge!(dict, Dict(f => string(getfield(mach.model, f))))
    end
    return Dict{Any,Any}(string(mach.model) => dict)
end

function write_results(outdir, results)
    run(`mkdir -p $outdir`)
    all_results = []
    for r in eachindex(results)
        for res in results[r]
            push!(all_results, res)
        end
    end
    full_results = vcat(all_results...)
    CSV.write(joinpath(outdir, "fit.csv"), full_results)
end

function batch(model_factory, model_name, reps=64)
    emb_file_dir = datadir("feature_embedding")
    embedding_files = readdir(emb_file_dir)
    
    outdir = datadir("model_fit", model_name)

    #run(`mkdir -p $outdir`)
    #run(`mkdir -p $(joinpath(outdir))`)

    lk = ReentrantLock()

    for (i, f) in enumerate(embedding_files)
        results = Any[[] for i in 1:reps]
        modelname = string(split(f,".")[1])

        Threads.@threads for r in 1:reps
            @info "Thread: $(Threads.threadid()), File: $f, Feature Set: $i/$(length(embedding_files))"
            model_instance = model_factory()
            df = CSV.read(joinpath(emb_file_dir, f), DataFrame)
            prediction_df, cv = fit_model(
                df, model_instance;
            )

            model_dir = joinpath(outdir, modelname, "predictions")
            run(`mkdir -p $model_dir`)

            lock(lk) do
                CSV.write(joinpath(model_dir, "$r.csv"), prediction_df)
                push!(results[Threads.threadid()], cv)    
            end 
        
        end


        write_results(
            joinpath(outdir, modelname), results
        )
    end
end
=#