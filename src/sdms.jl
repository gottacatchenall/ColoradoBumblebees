using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using Dates
using EvoTrees
using JSON
using CairoMakie

const SDT = SpeciesDistributionToolkit

# ============================================================================
# CONFIGURATION CONSTANTS
# ============================================================================

const EARTH_SYSTEM_MODELS = [
    GFDL_ESM4, 
    IPSL_CM6A_LR, 
    MPI_ESM1_2_HR, 
    MRI_ESM2_0, 
    UKESM1_0_LL
]

const CLIMATE_SCENARIOS = [
    SSP126,
    SSP370,
    SSP585 
]

const FUTURE_TIMESPANS = [
    Year(2011) => Year(2040),
    Year(2041) => Year(2070), 
    Year(2071) => Year(2100)
]


const BOUNDING_BOX = (left=-109.7, right=-101.8, bottom=34.5, top=42.5)

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

"""
    parse_occurrence_from_row(row)

Convert a CSV row into an Occurrence object with standardized fields.
"""
function parse_occurrence_from_row(row)
    OccurrencesInterface.Occurrence(;
        presence = row.occurrenceStatus == "PRESENT",
        what = row.species,
        when = DateTime(replace(row.eventDate, "Z" => "")),
        where = (row.decimalLongitude, row.decimalLatitude),
    )
end

"""
    load_occurrence_data(data_directory)

Load GBIF occurrence records and taxonomy data, returning a vector of Occurrence objects.
Joins species keys with species names from the taxonomy file.
"""
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

"""
    group_occurrences_by_species(occurrences, minimum_occurrences=50)

Split occurrence records into species-level groups, filtering to only include
species with at least `minimum_occurrences` records.
Returns a dictionary mapping species name to Occurrences object.
"""
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

"""
    load_baseline_climate_layers(chelsa_directory)

Load the 19 CHELSA bioclimatic variables for the baseline period.
Returns a vector of SDM layers.
"""
function load_baseline_climate_layers(chelsa_directory)
    baseline_directory = joinpath(chelsa_directory, "baseline")

    filenames = readdir(baseline_directory)
    
    chelsa_paths = [
        joinpath(baseline_directory, filenames[findfirst(fn -> occursin(pattern, fn), filenames)]) 
        for pattern in ["_bio$(i)_" for i in 1:19]
    ]
    
    # Load and convert to Float32
    layers = [SDMLayer(path; BOUNDING_BOX...) for path in chelsa_paths]
    return [Float32.(layer) for layer in layers]
end

"""
    load_future_climate_layers(base_directory, scenario, earth_system_model, timespan, bbox)

Load future climate projection layers for a specific scenario, ESM, and time period.
"""
function load_future_climate_layers(base_directory, scenario, earth_system_model, timespan, bbox)
    esm_string = replace(lowercase(string(earth_system_model)), "_" => "-")
    
    path = joinpath(
        base_directory,
        string(scenario),
        replace(string(earth_system_model), "_" => "-"),
    )
    
    start_year, end_year = timespan[1].value, timespan[2].value
    
    # Construct paths for all 19 bioclimatic variables
    layer_paths = [
        joinpath(
            path, 
            "chelsa_bio$(layer)_$(start_year)-$(end_year)_$(esm_string)_$(lowercase(string(scenario)))_v.2.1.tif"
        ) 
        for layer in 1:19
    ]
    
    return [Float32.(SDMLayer(path; bbox...)) for path in layer_paths]
end

# ============================================================================
# PSEUDOABSENCE GENERATION
# ============================================================================

"""
    generate_pseudoabsences(presence_layer, buffer_distance_km, class_balance_ratio)

Generate pseudoabsence points using a buffer around known presences.
- `buffer_distance_km`: minimum distance from presences in kilometers
- `class_balance_ratio`: ratio of absences to presences (1.0 = equal numbers)
"""
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

"""
    prepare_training_data(layers, presence_layer, absence_layer)

Combine environmental layers with presence/absence data into feature matrix X
and label vector y for model training.
"""
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

"""
    train_model(X_train, y_train)

Train a Gaussian process model using EvoTrees for species distribution modeling.
"""
function train_model(X_train, y_train)
    return EvoTrees.fit(
        EvoTreeGaussian(),
        x_train = X_train',
        y_train = y_train,
    )
end

"""
    predict_distribution(model, feature_matrix)

Generate predictions from a trained model.
"""
function predict_distribution(model, feature_matrix)
    return EvoTrees.predict(model, feature_matrix')
end

"""
    create_prediction_layer(model, environmental_layers)

Apply a trained model across the entire study area to create prediction and
uncertainty maps.
"""
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

"""
    calculate_evaluation_metrics(y_true, y_predicted, thresholds=0:0.001:1)

Calculate ROC-AUC, PR-AUC, TSS, and optimal threshold for binary classification performance.
"""
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
    
    # Find optimal threshold using True Skill Statistic
    optimal_threshold, threshold_index = findmax(trueskill.(confusion_matrices))
    
    return Dict(
        :prauc => pr_auc,
        :rocauc => roc_auc,
        :tss => trueskill(confusion_matrices[threshold_index]),
        :threshold => optimal_threshold
    )
end

"""
    aggregate_fold_statistics(fold_stats)

Aggregate statistics across k-fold cross-validation folds.
Returns mean and standard deviation for each metric.
"""
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

"""
    fit_baseline_sdm(
        occurrences, 
        environmental_layers; 
        pseudoabsence_buffer_distance=25.0, 
        class_balance=1.0, 
        k=4
    )

Fit a species distribution model using k-fold cross-validation.

# Arguments
- `occurrences`: Species occurrence records
- `environmental_layers`: Environmental predictor variables
- `pseudoabsence_buffer_distance`: Buffer around presences in km (default: 25)
- `class_balance`: Ratio of absences to presences (default: 1.0)
- `k`: Number of cross-validation folds (default: 4)

# Returns
Tuple of (models, range_map, uncertainty_map, statistics, presence_layer, absence_layer)
"""
function fit_baseline_sdm(
    occurrences,
    environmental_layers;
    pseudoabsence_buffer_distance = 25.0,
    class_balance = 1.0,
    k = 4
)
    # Create presence layer from occurrence points
    presence_layer = mask(environmental_layers[begin], occurrences)
    
    @info "    |-> Generating pseudoabsences..."
    absence_layer = generate_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)
    
    # Prepare training data
    features, labels = prepare_training_data(environmental_layers, presence_layer, absence_layer)
    
    # Create cross-validation folds
    fold_indices = SDeMo.kfold(labels, features, k=k)
    
    # Storage for results
    fold_statistics = []
    prediction_layers = []
    uncertainty_layers = []
    trained_models = []
    
    @info "    |-> Training $(k) cross-validation folds..."
    
    # Train and evaluate each fold
    for (train_idx, validation_idx) in fold_indices
        model = train_model(features[:, train_idx], labels[train_idx])
        
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

"""
    project_future_distribution(models, chelsa_directory, scenario, esm, timespan, bbox)

Project species distribution under future climate conditions using an ensemble of models.
"""
function project_future_distribution(models, chelsa_directory, scenario, esm, timespan, bbox)
    future_layers = load_future_climate_layers(chelsa_directory, scenario, esm, timespan, bbox)
    
    predictions = []
    uncertainties = []
    
    for model in models
        prediction, uncertainty = create_prediction_layer(model, future_layers)
        push!(predictions, prediction)
        push!(uncertainties, uncertainty)
    end
    
    return mean(predictions), mean(uncertainties)
end

"""
    create_esm_ensemble(models, chelsa_directory, scenario, timespan, bbox)

Create an ensemble prediction across multiple Earth System Models.
"""
function create_esm_ensemble(models, chelsa_directory, scenario, timespan, bbox)
    predictions = []
    uncertainties = []
    
    for (i_esm, esm) in enumerate(EARTH_SYSTEM_MODELS)
        @info "    |    |-> ESM: $esm [$(i_esm)/$(length(EARTH_SYSTEM_MODELS))]"
        
        prediction, uncertainty = project_future_distribution(
            models, chelsa_directory, scenario, esm, timespan, bbox
        )
        
        push!(predictions, prediction)
        push!(uncertainties, uncertainty)
    end
    
    return mean(predictions), mean(uncertainties)
end

# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

"""
    save_future_projections(output_directory, future_results)

Save future projection rasters to disk organized by scenario and time period.
"""
function save_future_projections(output_directory, future_results)
    for (scenario, scenario_dict) in future_results
        scenario_dir = joinpath(output_directory, string(scenario))
        
        for (timespan, timespan_dict) in scenario_dict
            @info timespan, timespan_dict
            # Format time period as string (e.g., "2011-2040")
            period_string = string(timespan[1].value) * "-" * string(timespan[2].value)
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

"""
    save_baseline_results(output_directory, baseline_results)

Save baseline SDM results including presence/absence data, predictions, and uncertainty.
"""
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

"""
    save_all_sdm_outputs(output_directory, species_name, results)

Save all SDM outputs including baseline, future projections, metrics, and visualization.
"""
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

"""
    create_species_distribution_models(
        data_directory, 
        output_directory, 
        chelsa_directory, 
        species_name; 
        k=5, 
        class_balance=1.0
    )

Complete workflow to create baseline and future species distribution models.

# Arguments
- `data_directory`: Path to occurrence data
- `output_directory`: Path for saving results
- `chelsa_directory`: Path to CHELSA climate data
- `species_name`: Name of target species
- `k`: Number of cross-validation folds (default: 5)
- `class_balance`: Ratio of pseudoabsences to presences (default: 1.0)
"""
function create_species_distribution_models(
    data_directory, 
    output_directory, 
    chelsa_directory,
    species_name;
    k = 5,
    class_balance = 1.0
)
    @info "Creating SDMs for $species_name"
    @info "="^70
    @info "\n"
    
    @info "[ 1/5 ] Loading baseline climate layers..."
    baseline_layers = load_baseline_climate_layers(chelsa_directory)
    
    @info "[ 2/5 ] Loading occurrence records..."
    all_occurrences = load_occurrence_data(data_directory)
    species_occurrences = group_occurrences_by_species(all_occurrences)
    target_occurrences = species_occurrences[species_name]
    
    @info "[ 3/5 ] Fitting baseline SDM..."
    models, range_map, uncertainty_map, statistics, presences, absences = fit_baseline_sdm(
        target_occurrences, 
        baseline_layers; 
        k = k,
        class_balance = class_balance
    )
    optimal_threshold = statistics[:threshold]["mean"]
    
    @info "[ 4/5 ] Projecting future distributions..."
    future_projections = Dict()
    
    for (i_scenario, scenario) in enumerate(CLIMATE_SCENARIOS)
        future_projections[scenario] = Dict()
        
        for (i_timespan, timespan) in enumerate(FUTURE_TIMESPANS)
            @info "    |$("-"^40)"
            @info "    |-> $scenario [$(i_scenario)/$(length(CLIMATE_SCENARIOS))]   |   $(timespan[1].value)-$(timespan[2].value) [$(i_timespan)/$(length(FUTURE_TIMESPANS))]"
            @info "    |$("-"^40)"
            # Create ensemble across ESMs
            prediction, uncertainty = create_esm_ensemble(
                models, chelsa_directory, scenario, timespan, BOUNDING_BOX
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