function fit_model(features, model)
    y, X, species_pairs = unpack(features, ==(:interaction), ∉([:bee, :plant]); rng=123)
    y = coerce(y, Multiclass{2})
    cv, pred = crossvalidation(model, X, y, species_pairs)
    
    thres = cv[!,:threshold][1]

    preddf = copy(features)
    preddf.prediction = pred
    preddf.thresholded_prediction .= pred .> thres

    preddf = preddf[!,[:bee,:plant,:prediction,:thresholded_prediction, :interaction]]
    return preddf,cv
end

function get_metadata(mach::Machine)
    fn = fieldnames(typeof(mach.model))
    dict = Dict()
    for f in fn
        dict = merge!(dict, Dict(f => string(getfield(mach.model, f))))
    end
    return Dict{Any,Any}(string(mach.model) => dict)
end

# brt = BRT(rng=MersenneTwister(Threads.threadid()))

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

    # NEW TODO: just drop anything with S2 in it, why waste time fitting when
    # you don't show it
    filter!(x-> !contains(x, "S2"), embedding_files)
    @info embedding_files



    # TODO comment after xgboost is done on cluster
    #complete_runs = readdir(outdir)
    #all_inputs = [string(split(s, ".")[1]) for s in embedding_files]
    #I = findall(s->s ∉ complete_runs,all_inputs)
    #embedding_files = embedding_files[I]
    #@info embedding_files

    #TODO uncomment after xgboost fix
    run(`mkdir -p $outdir`)
    run(`mkdir -p $(joinpath(outdir))`)

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
