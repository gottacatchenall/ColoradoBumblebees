function fit_and_project_sdms(; cluster=false)
    occurrence_df = load_occurrence_data()
    species = unique(occurrence_df.species)
    dfs = [filter(r->r.species == s, occurrence_df) for s in species]

    lk = ReentrantLock()

    Threads.@threads for i in eachindex(species)
        @info "Species: $(species[i])\tThread:$(Threads.threadid())"
        this_df = dfs[i]
        fit_sdm(species[i], this_df, lk, cluster=cluster)
        return
    end
end

function fit_sdm(species_name, occurrence_df, filelock; cluster=false)
    outdir = cluster ? joinpath("/scratch/mcatchen/BeeSDMs", species_name, "baseline") : joinpath(datadir("artifacts", "SDMs", species_name, "baseline"))
    mkpath(outdir)

    baseline_layers = load_chelsa_baseline() 

    occurrence_layer = SimpleSDMPredictor(zeros(Bool, size(baseline_layers[1])); SpeciesDistributionToolkit.boundingbox(baseline_layers[begin])... )
    convert_occurrence_to_tif!(occurrence_layer, occurrence_df)
    pres, abs, bgmask = get_pres_and_abs(occurrence_layer)

    X, y, pres_and_abs = get_features_and_labels(pres, abs, baseline_layers)
    Xtrain, Ytrain, Xtest, Ytest = test_train_split(X, y)

    model = fit_evotree(
        brt_params(GaussianBRT()); x_train=Xtrain, y_train=Ytrain, x_eval=Xtest, y_eval=Ytest
    )

    pred, uncert = predict_single_sdm(model, baseline_layers)
    fit_dict = compute_fit_stats_and_cutoff(pred, pres_and_abs, y)    

    lock(filelock) do
        SpeciesDistributionToolkit.save(joinpath(outdir, "prediction.tif"), pred)
        SpeciesDistributionToolkit.save(joinpath(outdir, "uncertainty.tif"), uncert)
        json_string = JSON.json(fit_dict)
        open(joinpath(outdir, "fit.json"), "w") do f
            JSON.print(f, json_string)
        end
    end 

    project_sdms(model, species_name, filelock, cluster=cluster)
end

function project_sdms(model, species, lk; cluster=false)
    ssps = [SSP1_26, SSP2_45, SSP3_70]
    years = TIMESPANS

    for ssp in ssps
        for year in years
            layers = load_layers(layer_paths)
            
            outdir = cluster ?  joinpath(datadir("/scratch/mcatchen/BeeSDMs", species, ssp, year)) : joinpath(datadir("artifacts", "SDMs", species, ssp, year))
            mkpath(outdir)

            pred, uncert = predict_single_sdm(model, layers)
            
            lock(lk) do 
                SpeciesDistributionToolkit.save(joinpath(outdir, "prediction.tif"), pred)
                SpeciesDistributionToolkit.save(joinpath(outdir, "uncertainty.tif"), uncert)
            end
        end
    end
end

function predict_single_sdm(model, layers)
    mat = zeros(Float32, length(layers), prod(size(layers[begin])))
    layers_to_matrix!(layers, mat)

    I = eachindex(layers[begin].grid)
    pred = EvoTrees.predict(model, mat')

    distribution = SimpleSDMPredictor(
        zeros(Float32, size(layers[begin])); SpeciesDistributionToolkit.boundingbox(layers[begin])...
    )
    distribution.grid[I] = pred[:, 1]

    uncertainty = SimpleSDMPredictor(zeros(Float32, size(layers[begin])); SpeciesDistributionToolkit.boundingbox(layers[begin])...)
    uncertainty.grid[I] = pred[:, 2]

    return rescale(distribution, (0, 1)), rescale(uncertainty, (0, 1))
end 