function make_sdms(species_name, occurrence_df; cluster=false, radius=400, buffer = 40, bias=2.)
    data = load_data()

    this_species_df = filter(x -> x.species == species_name, occurrence_df)

    baseline_layers = load_chelsa_baseline()
    occurrence_layer = SimpleSDMPredictor(zeros(Bool, size(baseline_layers[1])); SpeciesDistributionToolkit.boundingbox(baseline_layers[begin])...)
    convert_occurrence_to_tif!(occurrence_layer, this_species_df)

    pres, abs  = get_pres_and_abs(occurrence_layer, radius, buffer, bias)

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

    progbar = ProgressMeter.Progress(length(ssps) * length(years))


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

function get_pres_and_abs(presences, radius, buffer, bias)
    absences = background_thickening(presences, radius=radius, buffer=buffer, num_points=bias * sum(presences))
    replace!(absences, false => nothing)

    return presences, absences
end

function _check_bounds(template, i)
    sz = size(template)
    i[1] <= 0 && return false
    i[2] <= 0 && return false
    i[1] > sz[1] && return false
    i[2] > sz[2] && return false
    return true
end


function get_offset_mask(presences; buffer = 15, radius=100.)
    y, x = size(presences) # axes returned by size are flipped 
    bbox = SpeciesDistributionToolkit.boundingbox(presences)
    Δx = (bbox[:right] - bbox[:left]) / x   # how much a raster cell is in long
    Δy = (bbox[:top] - bbox[:bottom]) / y   # how much a raster cell is in lat
    centroid = [0.5(bbox[:right] + bbox[:left]), 0.5(bbox[:bottom] + bbox[:top])]
    
    _, radius_lat = SpeciesDistributionToolkit._known_point(centroid, radius, 0)
    radius_lon, _ = SpeciesDistributionToolkit._known_point(centroid, radius, π / 2)
    
    _, buffer_lat = SpeciesDistributionToolkit._known_point(centroid, buffer, 0)
    buffer_lon, _ = SpeciesDistributionToolkit._known_point(centroid, buffer, 0)

    max_offset_x, max_offset_y = Int32.(floor.([(radius_lon - centroid[1]) / Δx, (radius_lat - centroid[2]) / Δy]))
    radius_mask = OffsetArrays.OffsetArray(
        ones(Bool, 2max_offset_x + 1, 2max_offset_y + 1),
        -max_offset_x:max_offset_x,
        -max_offset_y:max_offset_y
    )
    for i in CartesianIndices(radius_mask)
        long_offset, lat_offset = abs.([i[1], i[2]]) .* [Δx, Δy]
        total_dist = sqrt(long_offset^2 + lat_offset^2)   # this is the dist in lat long

        dist_to_rad = min(radius_lon - centroid[1], radius_lat - centroid[2])
        dist_to_buffer = max(buffer_lon - centroid[1], buffer_lat - centroid[2])

        radius_mask[i] = total_dist <= dist_to_rad && total_dist > dist_to_buffer
    end
    offsets = CartesianIndices((-max_offset_x:max_offset_x, -max_offset_y:max_offset_y))

    radius_mask, offsets
end 

function _check_bounds(template, i)
    sz = size(template)
    i[1] <= 0 && return false
    i[2] <= 0 && return false
    i[1] > sz[1] && return false
    i[2] > sz[2] && return false
    return true
end

function get_thickening_layer(
    presences, 
    radius_mask,
    radius_offsets
)
    presence_only = mask(presences, presences)
    presence_idx = findall(x -> x == 1, presence_only.grid)
    background = similar(presences, Int32)
    background.grid[findall(!isnothing, background.grid)] .= 1
    for occ_idx in presence_idx
        within_radius_idx = radius_offsets .+ occ_idx
        for (i, cartesian_idx) in enumerate(within_radius_idx)
            if _check_bounds(background, cartesian_idx) && !isnothing(background.grid[cartesian_idx])
                if radius_mask[radius_offsets[i]] 
                    background.grid[cartesian_idx] += 1
                end
            end
        end
    end
    background ./ sum(background)
end

function background_thickening(presence; num_points = 400, radius = 150, buffer=30)
    radius_mask, radius_offsets = get_offset_mask(
        presence,
        radius = radius,
        buffer = buffer
    )
    l = get_thickening_layer(presence, radius_mask, radius_offsets)

    ps = [l.grid...]
    is = [CartesianIndices(l.grid)...]

    pas = [is[rand(Categorical(ps))] for _ in 1:num_points]

    background = similar(presence, Int32)
    background.grid[findall(!isnothing, background.grid)] .= 0

    background.grid[pas] .= 1

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
