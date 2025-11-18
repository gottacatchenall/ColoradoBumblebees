using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using Dates
using EvoTrees
using JSON
using CairoMakie
using Random


const SDT = SpeciesDistributionToolkit

# ============================================================================
# CONFIGURATION CONSTANTS
# ============================================================================

const EARTH_SYSTEM_MODELS = [
    "ACCESS-CM2",
    "BCC-CSM2-MR",
    "CMCC-ESM2",
    "EC-Earth3-Veg",
    "GISS-E2-1-G",
    "INM-CM5-0",
    "IPSL-CM6A-LR",
    "MIROC6",
    "MPI-ESM1-2-HR",
    "MRI-ESM2-0",
    "UKESM1-0-LL"
]

const CLIMATE_SCENARIOS = [
    "SSP126",
    "SSP245",
    "SSP370",
]

const FUTURE_TIMESPANS = [
    "2021-2040",
    "2041-2060",
    "2061-2080",
    "2081-2100"
]


const BOUNDING_BOX = (left=-109.7, right=-101.8, bottom=34.5, top=42.5)

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

function parse_occurrence_from_row(row)
    OccurrencesInterface.Occurrence(;
        presence = row.occurrenceStatus == "PRESENT",
        what = row.species,
        when = DateTime(replace(row.eventDate, "Z" => "")),
        where = (row.decimalLongitude, row.decimalLatitude),
    )
end

function load_occurrence_data(data_directory)
    # Load occurrence records
    gbif_data = CSV.read(joinpath(data_directory, "gbif.csv"), DataFrame)
    
    # Load taxonomy mapping
    taxa_data = CSV.read(joinpath(data_directory, "taxa.csv"), DataFrame)
    filter!(row->row.speciesKey ∈ taxa_data.species_key, gbif_data)

    species_key_to_name = Dict([row.species_key => row.species_name for row in eachrow(taxa_data)])
    
    # Map species keys to names
    gbif_data.species = [species_key_to_name[key] for key in gbif_data.speciesKey]
    
    return [parse_occurrence_from_row(row) for row in eachrow(gbif_data)]
end

function group_occurrences_by_species(occurrences, minimum_occurrences=50)
    # Extract genus and species only (ignore subspecies, variety, etc.)
    function extract_binomial_name(species_string)
        name_parts = split(species_string, " ")[1:2]
        return name_parts[1] * " " * name_parts[2]
    end
    
    unique_species = unique(map(x -> extract_binomial_name(x.what), occurrences))
    
    species_dict = Dict()
    for species_name in unique_species
        species_occurrences = Occurrences(
            filter(x -> extract_binomial_name(x.what) == species_name, SDT.elements(occurrences))
        )
        
        # Only include species with sufficient data
        if length(species_occurrences) > minimum_occurrences
            species_dict[species_name] = species_occurrences
        end
    end
    
    return species_dict
end

# ============================================================================
# ENVIRONMENTAL LAYER FUNCTIONS
# ============================================================================

function load_baseline_climate_layers(worldclim_dir)
    baseline_directory = joinpath(worldclim_dir, "baseline")

    filenames = readdir(baseline_directory)
    
    layer_paths = [
        joinpath(baseline_directory, filenames[findfirst(fn -> occursin(pattern, fn), filenames)]) 
        for pattern in ["_bio_$i.tif" for i in 1:19]
    ]
    
    # Load and convert to Float32
    layers = [SDMLayer(path; BOUNDING_BOX...) for path in layer_paths]
    return [Float32.(layer) for layer in layers]
end

function load_future_climate_layers(base_directory, scenario, earth_system_model, timespan)    
    path = joinpath(
        base_directory,
        scenario, 
        timespan
    )
        
    # Construct paths for all 19 bioclimatic variables
    layer_path =
        joinpath(
            path, 
            "wc2.1_30s_bioc_$(earth_system_model)_$(lowercase(scenario))_$timespan.tif"
        ) 
    
    layers =  [Float32.(SDMLayer(layer_path; bandnumber=i, BOUNDING_BOX...)) for i in 1:19]
    
    # Standardize to correct minor float errors in loading
    for l in layers
        l.x = (BOUNDING_BOX.left, BOUNDING_BOX.right)
        l.y = (BOUNDING_BOX.bottom, BOUNDING_BOX.top)
    end

    return layers
end

# ============================================================================
# PSEUDOABSENCE GENERATION
# ============================================================================

function generate_pseudoabsences(presence_layer, buffer_distance_km, class_balance_ratio)
    # Create background mask (all non-NA areas)
    background = pseudoabsencemask(DistanceToEvent, presence_layer)
    
    # Create buffer zone around presences (excluded from sampling)
    buffer = pseudoabsencemask(WithinRadius, presence_layer; distance = buffer_distance_km)
    
    # Apply no-data mask
    nodata_mask = nodata(buffer, true) 
    background.indices .= nodata_mask.indices
    
    # Sample pseudoabsences from background
    num_pseudoabsences = Int(round(class_balance_ratio * sum(presence_layer)))
    return backgroundpoints(background, num_pseudoabsences)
end

# ============================================================================
# MODEL TRAINING FUNCTIONS
# ============================================================================

function prepare_training_data(layers, presence_layer, absence_layer)
    # Extract environmental values at presence and absence locations
    X = Matrix(hcat([
        vcat(layer[findall(presence_layer)], layer[findall(absence_layer)]) 
        for layer in layers
    ]...)')
    
    # Create binary labels (1 = presence, 0 = absence)
    y = Bool.(vcat(
        [1 for _ in findall(presence_layer)], 
        [0 for _ in findall(absence_layer)]
    ))
    
    return X, y
end

function train_model(X_train, y_train; max_depth = 6)
    return EvoTrees.fit(
        EvoTreeGaussian(max_depth = max_depth),
        x_train = X_train',
        y_train = y_train,
    )
end

function predict_distribution(model, feature_matrix)
    return EvoTrees.predict(model, feature_matrix')
end

function create_prediction_layer(model, environmental_layers)
    # Initialize output layers
    prediction_layer = deepcopy(environmental_layers[begin])
    uncertainty_layer = deepcopy(environmental_layers[begin])
    
    # Extract features for all cells
    feature_matrix = Matrix(hcat([
        [layer[i] for layer in environmental_layers] 
        for i in eachindex(environmental_layers[1])
    ]...))
    
    # Generate predictions
    predictions = predict_distribution(model, feature_matrix)
    
    # Fill layers with predictions
    prediction_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 1]
    uncertainty_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 2]
    
    return prediction_layer, uncertainty_layer
end

# ============================================================================
# MODEL EVALUATION
# ============================================================================

function calculate_evaluation_metrics(y_true, y_predicted, thresholds=0:0.001:1)
    # Calculate confusion matrices across all thresholds
    confusion_matrices = [ConfusionMatrix(y_predicted .> t, y_true) for t in thresholds]
    false_positive_rates, true_positive_rates = fpr.(confusion_matrices), tpr.(confusion_matrices)
    
    # Calculate ROC-AUC using trapezoidal rule
    roc_dx = [reverse(false_positive_rates)[i] - reverse(false_positive_rates)[i-1] for i in 2:length(false_positive_rates)]
    roc_dy = [reverse(true_positive_rates)[i] + reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    roc_auc = sum(roc_dx .* (roc_dy ./ 2.0))
    
    # Calculate PR-AUC using trapezoidal rule
    precisions = ppv.(confusion_matrices)
    pr_dx = [reverse(true_positive_rates)[i] - reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    pr_dy = [reverse(precisions)[i] + reverse(precisions)[i-1] for i in 2:length(precisions)]
    pr_auc = sum(pr_dx .* (pr_dy ./ 2.0))
    
    # Find optimal threshold using MCC
    optimal_threshold, threshold_index = findmax(mcc.(confusion_matrices))
    
    return Dict(
        :prauc => pr_auc,
        :rocauc => roc_auc,
        :tss => trueskill(confusion_matrices[threshold_index]),
        :mcc => mcc(confusion_matrices[threshold_index]),
        :threshold => optimal_threshold
    )
end

function aggregate_fold_statistics(fold_stats)
    return Dict(
        metric => Dict("mean" => mean(values), "std" => std(values))
        for metric in keys(first(fold_stats)) 
        for values in [[fold[metric] for fold in fold_stats]]
    )
end

# ============================================================================
# BASELINE SDM FITTING
# ============================================================================

function fit_baseline_sdm(
    occurrences,
    environmental_layers;
    pseudoabsence_buffer_distance = 25.0,
    class_balance = 1.0,
    max_depth = 6,
    k = 4
)
    # Create presence layer from occurrence points
    presence_layer = mask(environmental_layers[begin], occurrences)
    
    @info "    |-> Generating pseudoabsences..."
    Random.seed!(123) # standardize across tuning tests
    absence_layer = generate_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)
    
    # Prepare training data
    features, labels = prepare_training_data(environmental_layers, presence_layer, absence_layer)
    
    # Create cross-validation folds
    Random.seed!(123) # standardize across tuning tests
    fold_indices = SDeMo.kfold(labels, features, k=k)
    
    # Storage for results
    fold_statistics = []
    prediction_layers = []
    uncertainty_layers = []
    trained_models = []
    
    @info "    |-> Training $(k) cross-validation folds..."
    
    # Train and evaluate each fold
    for (train_idx, validation_idx) in fold_indices
        model = train_model(features[:, train_idx], labels[train_idx]; max_depth = max_depth)
        
        # Evaluate on validation set
        validation_predictions = predict_distribution(model, features[:, validation_idx])[:, 1]
        push!(fold_statistics, calculate_evaluation_metrics(labels[validation_idx], validation_predictions))
        
        # Generate full prediction maps
        prediction, uncertainty = create_prediction_layer(model, environmental_layers)
        push!(prediction_layers, prediction)
        push!(uncertainty_layers, uncertainty)
        push!(trained_models, model)
    end
    
    # Aggregate results across folds
    aggregated_statistics = aggregate_fold_statistics(fold_statistics)
    optimal_threshold = aggregated_statistics[:threshold]["mean"]
    
    # Create binary range map using optimal threshold
    mean_prediction = mean(prediction_layers)
    range_map = Int.(mean_prediction .> optimal_threshold)
    uncertainty_map = mean(uncertainty_layers)
    
    return trained_models, range_map, uncertainty_map, aggregated_statistics, presence_layer, absence_layer
end

# ============================================================================
# FUTURE PROJECTION FUNCTIONS
# ============================================================================

function project_future_distribution(models, worldclim_directory, scenario, esm, timespan)
    future_layers = load_future_climate_layers(worldclim_directory, scenario, esm, timespan)
    
    predictions = []
    uncertainties = []
    
    for model in models
        prediction, uncertainty = create_prediction_layer(model, future_layers)
        push!(predictions, prediction)
        push!(uncertainties, uncertainty)
    end
    
    return mean(predictions), mean(uncertainties)
end

function create_esm_ensemble(models, worldclim_directory, scenario, timespan)
    predictions = []
    uncertainties = []
    
    for (i_esm, esm) in enumerate(EARTH_SYSTEM_MODELS)
        @info "    |    |-> ESM: $esm [$(i_esm)/$(length(EARTH_SYSTEM_MODELS))]"
        
        prediction, uncertainty = project_future_distribution(
            models, worldclim_directory, scenario, esm, timespan
        )
        
        push!(predictions, prediction)
        push!(uncertainties, uncertainty)
    end
    
    return mean(predictions), mean(uncertainties)
end

# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

function save_future_projections(output_directory, future_results)
    for (scenario, scenario_dict) in future_results
        scenario_dir = joinpath(output_directory, string(scenario))
        
        for (timespan, timespan_dict) in scenario_dict
            period_string = timespan
            period_dir = joinpath(scenario_dir, period_string)
            mkpath(period_dir)
            
            SDT.SimpleSDMLayers.save(
                joinpath(period_dir, "uncertainty.tif"),
                timespan_dict[:uncertainty]
            )
            SDT.SimpleSDMLayers.save(
                joinpath(period_dir, "range.tif"),
                timespan_dict[:range]
            )
        end
    end
end

function save_baseline_results(output_directory, baseline_results)
    mkpath(output_directory)
    SDT.SimpleSDMLayers.save(
        joinpath(output_directory, "presences.tif"),
        Float32.(baseline_results[:presences])
    )
    
    SDT.SimpleSDMLayers.save(
        joinpath(output_directory, "absences.tif"),
        Float32.(baseline_results[:absences])
    )
    SDT.SimpleSDMLayers.save(
        joinpath(output_directory, "uncertainty.tif"),
        baseline_results[:uncertainty]
    )
    SDT.SimpleSDMLayers.save(
        joinpath(output_directory, "range.tif"),
        baseline_results[:range]
    )
end

function save_all_sdm_outputs(output_directory, species_name, results)
    mkpath(output_directory)
    
    species_dir = joinpath(output_directory, species_name, "SDMs")
    mkpath(species_dir)
    
    # Save baseline 
    baseline_dir = joinpath(species_dir, "baseline")
    save_baseline_results(baseline_dir, results[:baseline])
    
    # Save futures
    futures_dir = joinpath(species_dir, "future")
    save_future_projections(futures_dir, results[:future])
    
    # Save evaluation metrics 
    open(joinpath(species_dir, "metrics.json"), "w") do file
        JSON.print(file, results[:metrics])
    end
    
    # Save baseline visualization
    fig = Figure()
    ax = Axis(fig[1, 1])
    heatmap!(ax, results[:baseline][:range])
    scatter!(ax, results[:baseline][:presences], color=:black, markersize=5)
    scatter!(ax, results[:baseline][:absences], color=:red, markersize=5)
    save(joinpath(species_dir, "baseline.png"), fig)

    return results
end

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

function create_species_distribution_models(
    data_directory, 
    output_directory, 
    worldclim_directory,
    species_name;
    k = 5,
    class_balance = 1.5,
    pseudoabsence_buffer_distance = 15.0,
)
    @info "Creating SDMs for $species_name"
    @info "="^70
    @info "\n"
    
    @info "[ 1/5 ] Loading baseline climate layers..."
    baseline_layers = load_baseline_climate_layers(worldclim_directory)
    
    @info "[ 2/5 ] Loading occurrence records..."
    all_occurrences = load_occurrence_data(data_directory)
    species_occurrences = group_occurrences_by_species(all_occurrences)
    target_occurrences = species_occurrences[species_name]
    
    @info "[ 3/5 ] Fitting baseline SDM..."
    models, range_map, uncertainty_map, statistics, presences, absences = fit_baseline_sdm(
        target_occurrences, 
        baseline_layers; 
        k = k,
        class_balance = class_balance,
        pseudoabsence_buffer_distance = pseudoabsence_buffer_distance,
    )
    optimal_threshold = statistics[:threshold]["mean"]
    
    @info "[ 4/5 ] Projecting future distributions..."
    future_projections = Dict()
    
    for (i_scenario, scenario) in enumerate(CLIMATE_SCENARIOS)
        future_projections[scenario] = Dict()
        
        for (i_timespan, timespan) in enumerate(FUTURE_TIMESPANS)
            @info "    |$("-"^40)"
            @info "    |-> $scenario [$(i_scenario)/$(length(CLIMATE_SCENARIOS))]   |   $(timespan) [$(i_timespan)/$(length(FUTURE_TIMESPANS))]"
            @info "    |$("-"^40)"
            # Create ensemble across ESMs
            prediction, uncertainty = create_esm_ensemble(
                models, worldclim_directory, scenario, timespan
            )
            
            future_projections[scenario][timespan] = Dict(
                :range => Int.(prediction .> optimal_threshold),
                :uncertainty => uncertainty
            )
        end
    end
    
    # Compile all results
    results = Dict(
        :baseline => Dict(
            :range => range_map,
            :uncertainty => uncertainty_map,
            :presences => presences,
            :absences => absences,
        ),
        :future => future_projections,
        :metrics => statistics,
    )
    
    @info "[ 5/5 ] Writing results..."
    save_all_sdm_outputs(output_directory, species_name, results)
end


function tune_hyperparameters(
    data_directory, 
    output_directory, 
    worldclim_directory,
    species_name;
    k = 5,
    class_balances = 0.5:0.5:3,
    pseudoabsence_buffer_distances = 5.0:5.0:25,
    max_depths = 4:2:10,
)
    
    @info "[ 1/3 ] Loading baseline climate layers..."
    baseline_layers = load_baseline_climate_layers(worldclim_directory)
    
    @info "[ 2/3 ] Loading occurrence records..."
    all_occurrences = load_occurrence_data(data_directory)
    species_occurrences = group_occurrences_by_species(all_occurrences)
    target_occurrences = species_occurrences[species_name]
    
    species_dir = joinpath(output_directory, species_name, "SDMs")
    mkpath(species_dir)

    results_df = DataFrame(
        class_balance = [],
        pseudoabsence_buffer_distance = [],
        max_depth = [],
        mcc = [],
        rocauc = []
    )

    cursor = 1
    num_treatments = prod(length.([class_balances, pseudoabsence_buffer_distances, max_depths]))

    @info "[ 3/3 ] Fitting models..."
    for class_balance in class_balances
        for pseudoabsence_buffer_distance in pseudoabsence_buffer_distances
            for max_depth in max_depths
                @info "    |-> Hyperparameter Set [$cursor / $num_treatments]..."
                _, _, _, statistics, _, _ = fit_baseline_sdm(
                    target_occurrences, 
                    baseline_layers; 
                    k = k,
                    class_balance = class_balance,
                    pseudoabsence_buffer_distance = pseudoabsence_buffer_distance,
                    max_depth = max_depth
                )

                push!(results_df, 
                    (
                        class_balance, 
                        pseudoabsence_buffer_distance, 
                        max_depth,
                        statistics[:mcc]["mean"],
                        statistics[:rocauc]["mean"]
                    )
                )

                # obnoxious way to check timing etc.
                CSV.write(joinpath(species_dir, "tuning.csv"), results_df)

                cursor += 1
            end 
        end
    end 

    CSV.write(joinpath(species_dir, "tuning.csv"), results_df)

    return results_df
end


function finish_tuning_hyperparameters(
    data_directory, 
    output_directory, 
    worldclim_directory,
    species_name;
    k = 5,
    class_balances = 0.5:0.5:3,
    pseudoabsence_buffer_distances = 5.0:5.0:25,
    max_depths = 4:2:10,
)
    
    @info "[ 1/3 ] Loading baseline climate layers..."
    baseline_layers = load_baseline_climate_layers(worldclim_directory)
    
    @info "[ 2/3 ] Loading occurrence records..."
    all_occurrences = load_occurrence_data(data_directory)
    species_occurrences = group_occurrences_by_species(all_occurrences)
    target_occurrences = species_occurrences[species_name]
    
    species_dir = joinpath(output_directory, species_name, "SDMs")
    mkpath(species_dir)



    results_df = CSV.read(joinpath(species_dir, "tuning.csv"), DataFrame)

    cursor = 1
    num_treatments = prod(length.([class_balances, pseudoabsence_buffer_distances, max_depths]))

    @info "[ 3/3 ] Fitting models..."
    for class_balance in class_balances
        for pseudoabsence_buffer_distance in pseudoabsence_buffer_distances
            for max_depth in max_depths

                if length(findall(
                    r-> r.class_balance == class_balance &&
                        r.pseudoabsence_buffer_distance == pseudoabsence_buffer_distance &&
                        r.max_depth == max_depth == max_depth,
                        eachrow(results_df)
                )) == 0 
                    @info "    |-> Hyperparameter Set [$cursor / $num_treatments]..."
                    _, _, _, statistics, _, _ = fit_baseline_sdm(
                        target_occurrences, 
                        baseline_layers; 
                        k = k,
                        class_balance = class_balance,
                        pseudoabsence_buffer_distance = pseudoabsence_buffer_distance,
                        max_depth = max_depth
                    )

                    push!(results_df, 
                        (
                            class_balance, 
                            pseudoabsence_buffer_distance, 
                            max_depth,
                            statistics[:mcc]["mean"],
                            statistics[:rocauc]["mean"]
                        )
                    )

                    CSV.write(joinpath(species_dir, "tuning.csv"), results_df)
                end
                cursor += 1
            end 
        end
    end 

    CSV.write(joinpath(species_dir, "tuning.csv"), results_df)

    return results_df
end