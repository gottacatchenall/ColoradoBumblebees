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


function _chelsa_dir_path(::Type{Timespan{S,E}}, s::Type{SC}) where {S,E,SC<:Scenario}
    startyear, endyear = S.value, E.value
    isbaseline = Timespan{S,E} == baseline()
    dirpath = isbaseline ? "baseline" : dirname(s)
    yrpath = isbaseline ? "" : "$startyear-$endyear"
    datadir("public", "chelsa", dirpath, yrpath)
end 

function load_chelsa(::Type{Timespan{S,E}}, s::Type{SC}) where {S,E, SC<:Scenario}
    dir = _chelsa_dir_path(Timespan{S,E}, s) 
    load_layers([joinpath(dir, x) for x in readdir(dir)])
end

function load_chelsa_baseline()
    load_chelsa(baseline(), Baseline)
end

load_chelsa(::Type{SC}, ::Type{Timespan{S,E}}) where {S,E, SC<:Scenario} = load_chelsa(Timespan{S,E}, SC)