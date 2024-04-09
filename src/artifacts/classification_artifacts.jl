function ColoradoBumblebees.save(bf::BatchFit; outdir=nothing)
    if isnothing(outdir)
        outdir = path(bf)
    end
    mkpath(outdir) 

    fit_stats_dir = joinpath(outdir, FIT_STATS_DIR) 
    mkpath(fit_stats_dir)

    prediction_dir = joinpath(outdir, PREDICTION_DF_DIR) 
    mkpath(prediction_dir)


    # write representation file if it doesn't exist
    # and add the pat hto metadata
    rep_file_path = _representation_file(bf)
    @info rep_file_path


    classification_model_metadata = _classification_model_to_dict(bf.fits[1].model)
    metadata = Dict(
        :classifier => classification_model_metadata,
        :representation_paths => rep_file_path, 
        :git_commit => gitdescribe(),
        :datetime => string(now()),
    )
    _write_json(joinpath(outdir, "metadata.json"), metadata)

    for (i,fit) in enumerate(bf.fits)
        _write_json(joinpath(fit_stats_dir, "fit$i.json"), fit.fit_stats)
        CSV.write(joinpath(prediction_dir, "fit$i.csv"), fit.predictions)
    end
end

function _load_classification_fit(path)
    metadata_path = joinpath(path, "metadata.json")    
    metadata = JSON.parsefile(metadata_path)

    classifier = _reconstruct_classifier(metadata["classifier"])
    @info metadata

    # Okay, perhaps this should be smarter about realizing where this was run
    # (e.g. local hostname or on cluster) and where it is currntly being run, to
    # make sure the path is valid.
    rep_path = metadata["representation_paths"]
    
    rep_path = [contains(x, "ColoradoBumblebees") ? x : joinpath("ColoradoBumblebees", x) for x in rep_path] 
    reps = ColoradoBumblebees.open.(rep_path)

    predictions = _load_predictions(path)
    fit_stats = _load_fit_stats(path)

    BatchFit([ClassificationFit(classifier, reps, predictions[i], fit_stats[i]) for i in eachindex(fit_stats)])
end

# What does this function do?
# Returns the paths to each species representation.
# If it is being used to write a batch fit (`bf`) and the representation hasn't
# been saved, it saves it.
function _representation_file(bf)
    function _create_if_nonexistant!(paths, r)
        if isnothing(r.path)
            @warn "This representation has not been saved. Saving representation before BatchFit."
            ColoradoBumblebees.save(r)
            push!(paths, path(r))
        else 
            push!(paths, r.path)
        end 
    end

    reps = bf.fits[begin].representation
    paths = []
    if typeof(reps) <: Vector
        for r in reps
            _create_if_nonexistant!(paths, r)
        end
    else
        _create_if_nonexistant!(paths, reps)
    end
    paths
end

function _classification_model_to_dict(model::ClassificationModel)
    params_dict = struct2dict(model)
    model_name = string(typeof(model))

    Dict(:model => model_name, :params => params_dict)
end


function _reconstruct_classifier(classifier_metadata)
    model_obj = _get_model_obj(classifier_metadata["model"])
    
    θ = Dict([Symbol(k) => v for (k,v) in classifier_metadata["params"]])
    model_obj(; θ...)
end



function _load_predictions(path)
    prediction_df_dir = joinpath(path, PREDICTION_DF_DIR)
    CSV.read.(joinpath.(prediction_df_dir, readdir(prediction_df_dir)), DataFrame)
end

function _load_fit_stats(path)
    fit_stats_dir = joinpath(path, FIT_STATS_DIR)
    _read_json.(joinpath.(fit_stats_dir, readdir(fit_stats_dir)))
end

_get_model_obj(model_string) = Dict("XGBoost" => XGBoost, "RandomForest" => RandomForest, "BoostedRegressionTree" => BoostedRegressionTree, "LogisticRegression"=>LogisticRegression, "Ensemble"=>Ensemble)[model_string]
