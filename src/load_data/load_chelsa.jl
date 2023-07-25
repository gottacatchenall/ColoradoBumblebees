function sort_pathlist(paths)
    I = sortperm(path_to_layer_num.(paths))
    paths[I]
end

function load_layer(path)
    df = CSV.read(path, DataFrame)
    bot,top = Float32.(extrema(df[!,:lat]))
    left, right = extrema(parse.(Float32, names(df)[2:end]))
    return SimpleSDMPredictor(Float32.(Matrix(df)), left=left, right=right, bottom=bot, top=top)
end

function load_layers(paths)
    layers = load_layer.(paths)
    for l in layers
        l.grid .= Float32.(l.grid)
    end
    layers
end

function load_baseline_layers()
    baseline_dir = joinpath(datadir("public", "chelsa", "baseline"))
    baseline_paths = sort_pathlist(joinpath.(baseline_dir, readdir(baseline_dir))) 
    baseline_layers = load_layers(baseline_paths)
end