# Todo 
# we want a function that loads occurrence data and loads the 19 decorrelated
# chelsa variables at that location for them 

function SpeciesDistributionToolkit.SimpleSDMLayers._match_longitude(layer::T, lon::K; lower::Bool=true) where {T <: SimpleSDMLayer, K <: AbstractFloat}
   # lon > 0. && layer.left <= lon <= layer.right || return nothing
   # lon < 0. && layer.left >= lon >= layer.right || return nothing

    lon == layer.left  && return 1
    lon == layer.right  && return size(layer, 2)
    relative = (lon - layer.left)/(layer.right - layer.left)
    fractional = relative * size(layer, 2)+1
    if lon in layer.left:2stride(layer,1):layer.right
        if abs(fractional - round(fractional)) < stride(layer, 1)
            fractional = round(fractional)
        end
        f = lower ? floor : ceil
        d = lower ? 1 : 0
        return min(f(Int64, fractional-d), size(layer, 2))
    else
        return min(floor(Int, fractional), size(layer, 2))
    end
end

function create_environmental_covariate_data()
    occurrence_df = load_occurrence_data()

    colprefix = "w_BIO" #"PCA_BIO"

    #current_layers = get_decorrelated_chelsa()
    current_layers = get_pca_chelsa()
    df = DataFrame(["species"=>[], ["$(colprefix)_$i" => [] for i in 1:length(current_layers)]...])
    for r in eachrow(occurrence_df)
        sp, lat, long = r
        push!(df.species, sp)
        for i in 1:length(current_layers)
            if !isnothing(current_layers[i][long,lat])
                push!(df[!, "$(colprefix)_$i"], current_layers[i][long,lat])
            else
                @info sp, long, lat
                push!(df[!, "$(colprefix)_$i"], NaN)
            end
        end
    end
    run(`mkdir -p $(datadir("public", "environment"))`)
    CSV.write(datadir("public", "environment", "covariates.csv"), df)
end

function get_decorrelated_chelsa()
    biolayers = ["BIO$i" for i in 1:19]
    current_layers = [convert(Float32,SimpleSDMPredictor(RasterData(CHELSA2, BioClim); layer=l, extent...)) for l in biolayers]
    I = common_Is(current_layers)
    mat = zeros(Float32,length(biolayers), length(I))
    w = fit_whitening(current_layers)
    current_decorrelated_layers = decorrelate_chelsa(current_layers, w, mat)
    standardize!(current_decorrelated_layers)
    current_decorrelated_layers
end

function get_pca_chelsa()
    biolayers = ["BIO$i" for i in 1:19]
    current_layers = [convert(Float32,SimpleSDMPredictor(RasterData(CHELSA2, BioClim); layer=l, extent...)) for l in biolayers]
    I = common_Is(current_layers)
    mat = zeros(Float32,length(biolayers), length(I))
    pca = fit_pca(current_layers)
    standardize!(pca_chelsa(current_layers, pca, mat))
end

function standardize!(current_layers)
    for l in current_layers
        μ, σ = mean(l.grid), std(l.grid)
        l.grid .= broadcast(x->(x-μ)/σ, l.grid)
    end
    current_layers
end

function common_Is(current_layers)
    Is = []
    for l in current_layers
        push!(Is, findall(x -> !isnothing(x) && !isnan(x), l.grid))
    end
    Is = unique(intersect(unique(Is)...))
end 

function fit_whitening(current_layers)
    Is = common_Is(current_layers)
    matrix = zeros(length(current_layers), length(Is))

    get_matrix_form!(current_layers, Is, matrix)
    matrix = convert.(Float32, matrix)

    @info "\t Fitting whitening..."
    w = MultivariateStats.fit(Whitening, matrix)
end

function fit_pca(current_layers)
    Is = common_Is(current_layers)
    matrix = zeros(length(current_layers), length(Is))

    get_matrix_form!(current_layers, Is, matrix)
    matrix = convert.(Float32, matrix)

    @info "\t Fitting PCA..."
    pca = MultivariateStats.fit(PCA, matrix)
end


function pca_chelsa(current_layers, pca, matrix)
    Is = common_Is(current_layers)
    get_matrix_form!(current_layers, Is, matrix)
    pca_matrix = MultivariateStats.transform(pca, matrix)
    new_layers = []
    for l in 1:size(pca_matrix)[1]
        tmp = convert(Float32,similar(current_layers[begin]))
        tmp.grid .= 0.
        #tmp.grid .= nothing
        tmp.grid[Is] .= pca_matrix[l,:]
        tmp
        push!(new_layers, tmp)
    end 
    new_layers
end

function decorrelate_chelsa(layers, w, matrix)
    Is = common_Is(layers)
    get_matrix_form!(layers, Is, matrix)
    decorrelated_matrix = MultivariateStats.transform(w, matrix)

    new_layers = []
    for l in 1:length(layers)
        tmp = convert(Float32,similar(layers[begin]))
        tmp.grid .= 0.
        #tmp.grid .= nothing
        tmp.grid[Is] .= decorrelated_matrix[l,:]
        tmp
        push!(new_layers, tmp)
    end 
    new_layers
end

function get_matrix_form!(layers, I, matrix)
    for l in 1:length(layers)
        for (ct,i) in enumerate(I)
            if isnothing(layers[l].grid[i]) || isnan(layers[l].grid[i])
                @info "fails, l:$l, i: $i, ct:$ct"
                return 
            else
                matrix[l,ct] = layers[l].grid[i]
            end
        end
    end 
    return matrix
end 


