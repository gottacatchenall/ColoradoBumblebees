function get_features_and_labels(presences, absences, climate_layers)
    presences = mask(presences, climate_layers[begin])
    absences = mask(absences, climate_layers[begin])
    coord_presence = keys(replace(presences, false => nothing))
    coord_absence = keys(replace(absences, false => nothing))
    coord = vcat(coord_presence, coord_absence)

    X = hcat([layer[coord] for layer in climate_layers]...)
    y = vcat(fill(1.0, length(coord_presence)), fill(0.0, length(coord_absence)))
    return X, y, coord
end

function convert_occurrence_to_tif!(occurrence_layer, occurrence_df)
    occurrence_layer.grid .= 0
    for r in eachrow(occurrence_df)
        long = SimpleSDMLayers._match_longitude(occurrence_layer, r.longitude)
        lat = SimpleSDMLayers._match_latitude(occurrence_layer, r.latitude)
    
        if !isnothing(long) && !isnothing(lat)
            i = CartesianIndex(long,lat)
            occurrence_layer.grid[i[2],i[1]] = 1.0   # long,lat is flipped in the grid representation

        end 
    end
end

function compute_fit_stats_and_cutoff(distribution, coords, y)
    cutoff = LinRange(extrema(distribution)..., 500)
    coords = convert(Vector{typeof(coords[begin])}, coords)
    idx = findall(!isnothing, coords)
    I = [SimpleSDMLayers._point_to_cartesian(distribution, c) for c in coords][idx]

    obs = y .> 0

    tp = zeros(Float64, length(cutoff))
    fp = zeros(Float64, length(cutoff))
    tn = zeros(Float64, length(cutoff))
    fn = zeros(Float64, length(cutoff))

    for (i, c) in enumerate(cutoff)
        prd = [distribution.grid[i] >= c for i in I]
        tp[i] = sum(prd .& obs)
        tn[i] = sum(.!(prd) .& (.!obs))
        fp[i] = sum(prd .& (.!obs))
        fn[i] = sum(.!(prd) .& obs)
    end

    tpr = tp ./ (tp .+ fn)
    fpr = fp ./ (fp .+ tn)
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0

    roc_dx = [reverse(fpr)[i] - reverse(fpr)[i - 1] for i in 2:length(fpr)]
    roc_dy = [reverse(tpr)[i] + reverse(tpr)[i - 1] for i in 2:length(tpr)]
    ROCAUC = sum(roc_dx .* (roc_dy ./ 2.0))

    thr_index = last(findmax(J))
    τ = cutoff[thr_index]

    return Dict(:rocauc => ROCAUC, :threshold => τ, :J => J[last(findmax(J))])
end

function layers_to_matrix!(climate_layers, mat)
    for (i, idx) in enumerate(eachindex(climate_layers[begin].grid))
        for l in eachindex(climate_layers)
            mat[l, i] = climate_layers[l].grid[idx]
        end
    end
end

function test_train_split(X, y, proportion=0.7)
    train_size = floor(Int, proportion * length(y))
    Itrain = StatsBase.sample(1:length(y), train_size; replace=false)
    Itest = setdiff(1:length(y), Itrain)
    Xtrain, Xtest = X[Itrain, :], X[Itest, :]
    Ytrain, Ytest = y[Itrain], y[Itest]
    return Xtrain, Ytrain, Xtest, Ytest
end