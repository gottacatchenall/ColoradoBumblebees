function clean_environmental_covariate_data()
    occurrence_df = load_occurrence_data()
    colprefix = "PCA_BIO"
    current_layers = get_pca_chelsa()
    df = DataFrame([
        "species" => [], ["$(colprefix)_$i" => [] for i in 1:length(current_layers)]...
    ])
    for r in eachrow(occurrence_df)
        _, sp, long, lat = r
        push!(df.species, sp)
        for i in eachindex(current_layers)
            push!(df[!, "$(colprefix)_$i"], current_layers[i][long, lat])
        end
    end
    filter!(r -> !any([isnothing(x) for x in r]), df)
    run(`mkdir -p $(datadir("public", "environment"))`)
    return CSV.write(datadir("public", "environment", "covariates.csv"), df)
end

function get_decorrelated_chelsa()
    current_layers = load_baseline_layers()
    I = common_Is(current_layers)
    mat = zeros(Float32, length(biolayers), length(I))
    w = fit_whitening(current_layers)
    current_decorrelated_layers = decorrelate_chelsa(current_layers, w, mat)
    standardize!(current_decorrelated_layers)
    return current_decorrelated_layers
end

function path_to_layer_num(path)
    s = split(path, "bio")
    parse(Int32, split(split(s[end], ".")[begin], "/")[end])
end

function get_pca_chelsa()
    current_layers = load_baseline_layers()
    I = common_Is(current_layers)
    mat = zeros(Float32, length(current_layers), length(I))
    pca = fit_pca(current_layers)
    return standardize!(pca_chelsa(current_layers, pca, mat))
end

function standardize!(current_layers)
    for l in current_layers
        μ, σ = mean(l.grid), std(l.grid)
        l.grid .= broadcast(x -> (x - μ) / σ, l.grid)
    end
    return current_layers
end

function common_Is(current_layers)
    Is = []
    for l in current_layers
        push!(Is, findall(x -> !isnothing(x) && !isnan(x), l.grid))
    end
    return Is = unique(intersect(unique(Is)...))
end

function fit_whitening(current_layers)
    Is = common_Is(current_layers)
    matrix = zeros(length(current_layers), length(Is))

    get_matrix_form!(current_layers, Is, matrix)
    matrix = convert.(Float32, matrix)

    @info "\t Fitting whitening..."
    return w = MultivariateStats.fit(Whitening, matrix)
end

function fit_pca(current_layers)
    Is = common_Is(current_layers)
    matrix = zeros(length(current_layers), length(Is))

    get_matrix_form!(current_layers, Is, matrix)
    matrix = convert.(Float32, matrix)

    @info "\t Fitting PCA..."
    return pca = MultivariateStats.fit(PCA, matrix)
end

function pca_chelsa(current_layers, pca, matrix)
    Is = common_Is(current_layers)
    get_matrix_form!(current_layers, Is, matrix)
    pca_matrix = MultivariateStats.transform(pca, matrix)
    new_layers = []
    for l in 1:size(pca_matrix)[1]
        tmp = convert(Float32, similar(current_layers[begin]))
        tmp.grid .= 0.0
        #tmp.grid .= nothing
        tmp.grid[Is] .= pca_matrix[l, :]
        tmp
        push!(new_layers, tmp)
    end
    return new_layers
end

function decorrelate_chelsa(layers, w, matrix)
    Is = common_Is(layers)
    get_matrix_form!(layers, Is, matrix)
    decorrelated_matrix = MultivariateStats.transform(w, matrix)

    new_layers = []
    for l in 1:length(layers)
        tmp = convert(Float32, similar(layers[begin]))
        tmp.grid .= 0.0
        #tmp.grid .= nothing
        tmp.grid[Is] .= decorrelated_matrix[l, :]
        tmp
        push!(new_layers, tmp)
    end
    return new_layers
end

function get_matrix_form!(layers, I, matrix)
    for l in 1:length(layers)
        for (ct, i) in enumerate(I)
            if isnothing(layers[l].grid[i]) || isnan(layers[l].grid[i])
                @info "fails, l:$l, i: $i, ct:$ct"
                return nothing
            else
                matrix[l, ct] = layers[l].grid[i]
            end
        end
    end
    return matrix
end
