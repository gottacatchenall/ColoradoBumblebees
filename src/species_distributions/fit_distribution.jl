function make_sdms(species_name, occurrence_df; cluster=false, pa_buffer_distance=20.)
    data = load_data()
    
    baseline_layers = load_chelsa_baseline() 
    occurrence_layer = SimpleSDMPredictor(zeros(Bool, size(baseline_layers[1])); SpeciesDistributionToolkit.boundingbox(baseline_layers[begin])... )
    convert_occurrence_to_tif!(occurrence_layer, occurrence_df)

    pres, abs, bgmask = get_pres_and_abs(occurrence_layer; pa_buffer_distance=pa_buffer_distance)

    X, y, pres_and_abs = get_features_and_labels(pres, abs, baseline_layers)
    Xtrain, Ytrain, Xtest, Ytest = test_train_split(X, y)

    model = fit_evotree(
        brt_params(BoostedRegressionSDM()); 
        x_train=Xtrain, 
        y_train=Ytrain, 
        x_eval=Xtest, 
        y_eval=Ytest
    )

    prob, uncert = predict_single_sdm(model, baseline_layers)
    fit_dict = compute_fit_stats_and_cutoff(prob, pres_and_abs, y)    
    
    species = contains(species_name, "Bombus") ? bee(data, species_name) : plant(data, species_name)

    baseline_sdm = SpeciesDistribution(species, prob, uncert, fit_dict, baseline(), Baseline)

    future_projections = project_sdms(model, species, fit_dict, cluster=cluster)

    [baseline_sdm, future_projections...]
end

function project_sdms(model, species, fit_dict; cluster=false)
    ssps = [SSP1_26, SSP2_45, SSP3_70]
    years = TIMESPANS[2:end]

    distributions = []

    progbar = ProgressMeter.Progress(length(ssps)*length(years))
          

    for ssp in ssps
        for year in years
            layers = load_chelsa(year, ssp)
            prob, uncert = predict_single_sdm(model, layers)
            

            push!(distributions, SpeciesDistribution(species, prob, uncert, fit_dict, year, ssp))
            ProgressMeter.next!(
                progbar; showvalues=[(Symbol("SSP"), ssp), (Symbol("Years"), year)]
            )
        end
    end
    distributions
end

function get_pres_and_abs(presences; pa_buffer_distance = 10.)
    @time buffer = _new_pa(presences; distance = pa_buffer_distance)
    bgmask = (.!buffer)
    absences = SpeciesDistributionToolkit.sample(bgmask, floor(Int, 0.5sum(presences)))
    replace!(absences, false => nothing)

    return presences, absences, bgmask
end

function _check_bounds(template, i)
    sz = size(template)
    i[1] <= 0 && return false
    i[2] <= 0 && return false
    i[1] > sz[1] && return false
    i[2] > sz[2] && return false
    return true
end

function _new_pa(
    presences::T;
    distance::Number = 100.0,
) where {T <: SimpleSDMLayer}
    presence_only = mask(presences, presences)
    presence_idx = findall(x->x==1,presence_only.grid)

    background = similar(presences, Bool)
    background.grid .= true

    y,x = size(presences)
    bbox = boundingbox(presences)
    Δx = (bbox[:right] - bbox[:left])/ x   # how much a given cell is in long
    Δy = (bbox[:top] - bbox[:bottom])/ y   # how much a given cell is in lat
    
    lon = zeros(Float64, 2)
    lat = zeros(Float64, 2)
    for (i, angl) in enumerate((0:1) / 4)
        α = deg2rad(360.0angl)
        lon[i], lat[i] = SpeciesDistributionToolkit._known_point([0.0, 0.0], distance, α)
    end

    max_cells_x, max_cells_y =  Int32.(floor.([lon[2] / Δx, lat[1] / Δy])) 

    radius_mask = OffsetArray(ones(Bool, 2max_cells_x+1, 2max_cells_y+1), -max_cells_x:max_cells_x, -max_cells_y:max_cells_y) 

    for i in CartesianIndices(radius_mask)
        long_offset, lat_offset = abs.([i[1], i[2]]) .* [Δx, Δy]

        total_dist = sqrt(long_offset^2 + lat_offset^2)
        radius_mask[i] = total_dist <= max(lon[2], lat[1])
    end

    # Mask radius around each presence point 
    for occ_idx in presence_idx
        offs = CartesianIndices((-max_cells_x:max_cells_x, -max_cells_y:max_cells_y))
        I = offs .+ occ_idx

        for (i, idx) in enumerate(I)
            if _check_bounds(background, idx) && radius_mask[offs[i]]
                background[idx] = false
            end 
        end 
    end
    
    return background
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