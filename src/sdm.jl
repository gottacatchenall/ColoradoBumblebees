using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using Dates
using EvoTrees
using JSON
using CairoMakie

const SDT = SpeciesDistributionToolkit

function occurrence_from_row(row)
    OccurrencesInterface.Occurrence(;
        presence = row.occurrenceStatus == "PRESENT",
        what = row.species,
        when = DateTime(replace(row.eventDate, "Z" => "")),
        where = (row.decimalLongitude, row.decimalLatitude),
    )
end

function get_occurrences(data_dir)
    gbif_df = CSV.read(joinpath(data_dir, "gbif.csv"), DataFrame)
    taxa_df = CSV.read(joinpath(data_dir, "taxa.csv"), DataFrame)
    key2name = Dict([x.species_key => x.species_name for x in eachrow(taxa_df)])
    gbif_df.species = [key2name[x] for x in gbif_df.speciesKey]
    return [occurrence_from_row(r) for r in eachrow(gbif_df)]
end


function write_sdm_artifacts(artifact_dir, species_name, results)
    mkpath(artifact_dir)

    output_dir = joinpath(artifact_dir, species_name)
    mkpath(output_dir)

    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "presences.tif"),
        Float32.(results[:presences])
    )      
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "absences.tif"),
        Float32.(results[:absences])
    )      
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "uncertainty.tif"),
        results[:uncertainty]
    )        
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "prediction.tif"),
        results[:prediction]
    )        
    open(joinpath(output_dir, "metrics.json"), "w") do f
        JSON.print(f, results[:metrics])
    end

    f = Figure()
    ax = Axis(f[1,1])
    heatmap!(ax, results[:prediction])
    scatter!(ax, results[:presences], color=:black, markersize=5)
    scatter!(ax, results[:absences], color=:red, markersize=5)
    save(joinpath(output_dir, "plot.png"), f)
end


function split_occurrences_into_species(occs, min_occs=50)
    function _transform_to_genus_species(str)
        x = split(str, " ")[1:2]
        return x[1]*" "*x[2]
    end
    unique_species = unique(map(
        x->_transform_to_genus_species(x.what),
        occs)
    )

    dict = Dict()
    for sp in unique_species
        these_occs = Occurrences(filter(x->_transform_to_genus_species(x.what) == sp, SDT.elements(occs)))
        if length(these_occs) > min_occs
            dict[sp] = these_occs
        end
    end
    return dict
end

function get_pseudoabsences(
    presence_layer, 
    buffer_distance,
    class_balance
)
    background = pseudoabsencemask(DistanceToEvent, presence_layer)
    buffer = pseudoabsencemask(WithinRadius, presence_layer; distance = buffer_distance)
    ndmask = nodata(buffer, true) 
    background.indices .= ndmask.indices

    pseudoabs_layer = backgroundpoints(background, Int(round(class_balance*sum(presence_layer))))
    return pseudoabs_layer
end

function get_features_and_labels(layers, presence_layer, absence_layer)
    X = Matrix(hcat([vcat(l[findall(presence_layer)], l[findall(absence_layer)]) for l in layers]...)')
    y = Bool.(vcat([1 for _ in findall(presence_layer)], [0 for _ in findall(absence_layer)]))
    X, y
end

function predict_sdm(model, X)
    EvoTrees.predict(model, X')
end

function fit_fold(X, y)
    model = EvoTrees.fit_evotree(
        EvoTreeGaussian(),
        x_train=X',
        y_train=y,
    )
end

function get_fit_stats(y_true, y_predict, τ = 0:0.001:1)
    cms = [ConfusionMatrix(y_predict .> t, y_true) for t in τ]
    FPRs, TPRs = fpr.(cms), tpr.(cms)

    roc_dx = [reverse(FPRs)[i] - reverse(FPRs)[i - 1] for i in 2:length(FPRs)]
    roc_dy = [reverse(TPRs)[i] + reverse(TPRs)[i - 1] for i in 2:length(TPRs)]
    ROCAUC = sum(roc_dx .* (roc_dy ./ 2.0))

    precisions = ppv.(cms)
    pr_dx = [reverse(TPRs)[i] - reverse(TPRs)[i - 1] for i in 2:length(TPRs)]
    pr_dy = [reverse(precisions)[i] + reverse(precisions)[i - 1] for i in 2:length(precisions)]
    PRAUC = sum(pr_dx .* (pr_dy ./ 2.0))

    threshold, threshold_idx = findmax(trueskill.(cms))

    return Dict(
        :prauc => PRAUC,
        :rocauc => ROCAUC,
        :tss => trueskill(cms[threshold_idx]),
        :threshold => threshold
    )
end

function get_prediction_layer(model, layers)
    prediction, uncertainty = deepcopy(layers[begin]),  deepcopy(layers[begin])
    prediction_mat = predict_sdm(model, Matrix(hcat([[l[i] for l in layers] for i in eachindex(layers[1])]...)))
    
    prediction.grid[findall(prediction.indices)] .= prediction_mat[:,1]
    uncertainty.grid[findall(prediction.indices)] .= prediction_mat[:,2]
    return prediction, uncertainty
end 

function aggregate_dicts(dicts)
    Dict(k => Dict("mean" => mean(values), "std" => std(values))
        for k in keys(first(dicts)) 
        for values in [[d[k] for d in dicts]]
    )
end


function fit_sdm(
    occurrences,
    layers;
    pseudoabsence_buffer_distance = 25.,
    class_balance = 1., # high val -> more negatives
    k = 4
)
    @info "Starting SDM for $(occurrences[1].what)"
    
    presence_layer = mask(layers[begin], occurrences)
    absence_layer = get_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)
    features, labels = get_features_and_labels(layers, presence_layer, absence_layer)
    fold_indices = SDeMo.kfold(labels, features, k=k)

    fit_stats = []
    prediction_layers = []
    uncertainty_layers = []

    # Fit folds
    for (train_idx, val_idx) in fold_indices
        model = fit_fold(features[:, train_idx], labels[train_idx])
        y_predict = predict_sdm(model, features[:, val_idx])[:, 1] # first col is prediction
        push!(fit_stats, get_fit_stats(labels[val_idx], y_predict))

        p, u = get_prediction_layer(model, layers)
        push!(prediction_layers, p)
        push!(uncertainty_layers, u)
    end

    return mean(prediction_layers), mean(uncertainty_layers), aggregate_dicts(fit_stats), presence_layer, absence_layer
end 



