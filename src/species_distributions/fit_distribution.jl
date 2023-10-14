function make_sdms(species_name, occurrence_df; cluster=false, pa_buffer_distance=15.)
    data = load_data()

    this_species_df = filter(x->x.species == species_name, occurrence_df)
    
    baseline_layers = load_chelsa_baseline() 
    occurrence_layer = SimpleSDMPredictor(zeros(Bool, size(baseline_layers[1])); SpeciesDistributionToolkit.boundingbox(baseline_layers[begin])... )
    convert_occurrence_to_tif!(occurrence_layer, this_species_df)

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

function get_pres_and_abs(presences; radius=50., pa_buffer_distance = 3.)
    valid_regions = pseudoabsencemask(WithinRadius, presences; distance = radius)
    
    nothing_idx = findall(isnothing, valid_regions.grid)
    valid_regions.grid[nothing_idx] .= false
    
    too_close = _new_pa(presences; distance = pa_buffer_distance)

    bgmask = .!valid_regions .& too_close
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
        
    function _check_bounds(template, i)
        sz = size(template)
        i[1] <= 0 && return false
        i[2] <= 0 && return false
        i[1] > sz[1] && return false
        i[2] > sz[2] && return false
        return true
    end

    presence_only = mask(presences, presences)
    presence_idx = findall(x->x==1,presence_only.grid)


    y, x = size(presences) # axes returned by size are flipped 
    bbox = boundingbox(presences)
    Δx = (bbox[:right] - bbox[:left])/ x   # how much a raster cell is in long
    Δy = (bbox[:top] - bbox[:bottom])/ y   # how much a raster cell is in lat

    # It's reasonable to use the centroid of the raster as a basis and use the
    # same sliding mask for each point.
    centroid = [0.5(bbox[:right]+bbox[:left]), 0.5(bbox[:bottom]+bbox[:top])]

    # However, if the extent is large enough, this assumption might break down
    # and the size of the sliding window should be calculated for each
    # occurrence (or for some subbdivision of the raster into subsections).

    # The slowest (but most accurrate) version would create an offset mask for
    # each occurrence. This is still likely faster than the current version
    # method of filtering coordinates, but may be significantly slower than the
    # simpler methods.  

    _, lat = SpeciesDistributionToolkit._known_point(centroid, distance, 0)
    lon, _ = SpeciesDistributionToolkit._known_point(centroid, distance, π/2)

    # magnitude of the maximum offset in lon/lat from raster centroid in each direction 
    max_offset_x, max_offset_y =  Int32.(floor.([(lon-centroid[1]) / Δx, (lat-centroid[2]) / Δy])) 

    radius_mask = OffsetArrays.OffsetArray(
        ones(Bool, 2max_offset_x+1, 2max_offset_y+1), 
        -max_offset_x:max_offset_x, 
        -max_offset_y:max_offset_y
    ) 
    for i in CartesianIndices(radius_mask)
        long_offset, lat_offset = abs.([i[1], i[2]]) .* [Δx, Δy]
        total_dist = sqrt(long_offset^2 + lat_offset^2)
        radius_mask[i] = total_dist <= min(lon-centroid[1], lat-centroid[2])  
    end

    # Mask radius around each presence point 
    offsets = CartesianIndices((-max_offset_x:max_offset_x, -max_offset_y:max_offset_y))

    # 3 states 
    # State 1: unseen
    # State 2: seen, still a candidate bg point
    # State 3: seen, no longer a cnadidate bg point
    UNSEEN, CANDIDATE, UNAVAILABLE = 1, 2, 3 

    background = similar(presences, Int32)
    background.grid[findall(!isnothing, background.grid)] .= UNSEEN

    for occ_idx in presence_idx
        within_radius_idx = offsets .+ occ_idx

        for (i, cartesian_idx) in enumerate(within_radius_idx)
            if _check_bounds(background, cartesian_idx) && !isnothing(background.grid[cartesian_idx])
                if background.grid[cartesian_idx] != UNAVAILABLE && !radius_mask[offsets[i]] 
                    background.grid[cartesian_idx] = CANDIDATE
                else
                    background.grid[cartesian_idx] = UNAVAILABLE
                end 
            end
        end 
    end

    background.grid[findall(isequal(CANDIDATE), background.grid)] .= true
    background.grid[findall(isequal(UNAVAILABLE), background.grid)] .= false

    return convert(Bool, background)
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