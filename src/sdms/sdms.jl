#smolext = (bottom=35.0, top=40.0, left=-108.0, right=-105.0)
#using CSV, DataFrames, MultivariateStats, EvoTrees, SpeciesDistributionToolkit, JSON

Base.@kwdef struct GaussianBRT
    loss = :gaussian
    metric = :gaussian
    nrounds = 100
    nbins = 100
    λ = 0.0
    γ = 0.0
    η = 0.1
    max_depth = 7
    min_weight = 1.0
    rowsample = 0.5
    colsample = 1.0
end

function brt_params(gbrt::GaussianBRT)
    return EvoTreeGaussian(; (v => getfield(gbrt, v) for v in fieldnames(typeof(gbrt)))...)
end

function make_sdms(data, gbrt::GaussianBRT; cluster=false)
    !cluster && @warn "Running on small extent"
    extent = cluster ? EXTENT : (bottom=35.0, top=40.0, left=-108.0, right=-105.0)
    current_layers = [
        convert(
            Float32, SimpleSDMPredictor(RasterData(CHELSA2, BioClim); layer=l, extent...)
        ) for l in BIOLAYERS
    ]

    mat = zeros(Float32, length(current_layers), prod(size(current_layers[begin])))
    w = fit_whitening(current_layers, mat)
    current_layers = decorrelate_chelsa(current_layers, w, mat)

    # TESTING
    #@warn "Still using 6 species"
    #species = plants(data)[1:6]
    @warn "Running on all species"
    species = vcat(bees(data)..., plants(data)...)

    occurrence_layer = convert(Bool, similar(current_layers[begin]))
    occurrence_layer.grid .= 0

    scen = scenarios()
    !cluster && @warn "Still only running on ssp370"
    baseline, futures = scen[1], cluster ? scen[2:end] : scen[5:7]

    models = fit_all_sdms(
        species, current_layers, occurrence_layer, gbrt, baseline, mat; cluster=cluster
    )
    return project_future_sdms(species, models, current_layers, w, futures, mat, extent; cluster=cluster)
end

function project_future_sdms(species, models, current_layers, w, futures, mat, extent; cluster=false)
    for f in futures
        current_layers = load_climate_layers(f, extent)
        current_layers = decorrelate_chelsa(current_layers, w, mat)
        project_sdms_single_future(f, species, models, current_layers, mat; cluster=cluster)
    end
end

function project_sdms_single_future(scenario, species, models, current_layers, mat; cluster=cluster)
    sdm_out_paths = [get_sdm_path(scenario, sp; cluster=cluster) for sp in species]
    uncertainty_out_paths = [get_uncertainty_path(scenario, sp; cluster=cluster) for sp in species]
    sdm_dirs = [get_sdm_dir(scenario, sp; cluster=cluster) for sp in species]

    sdms, uncertainties = [], []

    Threads.@threads for THREAD_ID in 1:Threads.nthreads()
        sp = species[THREAD_ID]
        outdir = sdm_dirs[THREAD_ID]
        run(`mkdir -p $outdir`)
        thismodel = models[sp]
        predicted_sdm, predicted_uncertainty = predict_sdm(current_layers, thismodel, zeros(size(mat)))
        push!(sdms, predicted_sdm)
        push!(uncertainties, predicted_uncertainty)
    end

    for i in eachindex(sdms)
        SpeciesDistributionToolkit.save(sdm_out_paths[i], sdms[i]; driver="GTiff")
        SpeciesDistributionToolkit.save(
            uncertainty_out_paths[i], uncertainties[i]; driver="GTiff"
        )
    end
end

function fit_all_sdms(
    species, climate_layers, template_layer, gbrt, baseline, mat; cluster=false
)
    models = Dict()

    sdm_out_paths = [get_sdm_path(baseline, sp; cluster=cluster) for sp in species]
    uncertainty_out_paths = [get_uncertainty_path(baseline, sp; cluster=cluster) for sp in species]
    sdm_dirs = [get_sdm_dir(baseline, sp; cluster=cluster) for sp in species]

    sdms = []
    uncertainties = []
    fitstats = []

    for THREAD_ID in eachindex(species)
        sp = species[THREAD_ID]
        outdir = sdm_dirs[THREAD_ID]
        @info sp
        run(`mkdir -p $outdir`)
       
        #occurrence_layer = similar(template_layer)
        #occurrence_layer.grid .= 0
        #convert_occurrence_to_tif!(sp, occurrence_layer)
        convert_occurrence_to_tif!(sp, template_layer)
        model, coords, labels = fit_sdm(climate_layers, template_layer, gbrt)
        
        # THIS ISNT MEMORY SAFE w/ THREADS FIND A NEW WAY
        merge!(models, Dict(sp => model))

        #predicted_sdm, predicted_uncertainty = predict_sdm(climate_layers, model, zeros(size(mat)))
        predicted_sdm, predicted_uncertainty = predict_sdm(climate_layers, model, mat)

        write_stats(
            compute_fit_stats_and_cutoff(predicted_sdm, coords, labels), joinpath(sdm_dirs[THREAD_ID], "fit.json")
        )
        SpeciesDistributionToolkit.save(
            sdm_out_paths[THREAD_ID], predicted_sdm; driver="GTiff"
        )
        SpeciesDistributionToolkit.save(
            uncertainty_out_paths[THREAD_ID], predicted_uncertainty; driver="GTiff"
        )
        #push!(sdms, predicted_sdm)
        #push!(uncertainties, predicted_uncertainty)
       #println("About to compute stats for $sp with coords: $(typeof(coords)), $(length(coords))")
       # push!(fitstats, compute_fit_stats_and_cutoff(predicted_sdm, coords, labels))
    end

    #= for i in eachindex(sdms)
        write_stats(
            fitstats[i], joinpath(sdm_dirs[i], "fit.json")
        )
        SpeciesDistributionToolkit.save(
            sdm_out_paths[i], sdms[i]; driver="GTiff"
        )
        SpeciesDistributionToolkit.save(
            uncertainty_out_paths[i], uncertainties[i]; driver="GTiff"
        )
    end =#
    return models
end

_scenario_to_projection(scenario, model=GFDL_ESM4) = Projection(scenario.ssp, model)
function _scenario_to_year_pair(scenario)
    return scenario.years.startyear => scenario.years.endyear
end

function load_climate_layers(scenario, extent)
    return tmp = [
        SimpleSDMPredictor(
            RasterData(CHELSA2, BioClim),
            _scenario_to_projection(scenario);
            timespan=_scenario_to_year_pair(scenario),
            layer=l,
            extent...,
        ) for l in BIOLAYERS
    ]
end

function SimpleSDMDatasets.source(data::RasterData{CHELSA2,BioClim}; layer="BIO1")
    var_code = (layer isa Integer) ? layer : findfirst(isequal(layer), layers(data))
    root = "https://envicloud.wsl.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/"
    stem = "CHELSA_bio$(var_code)_1981-2010_V.2.1.tif"
    return (url=root * stem, filename=lowercase(stem), outdir=destination(data))
end

function SimpleSDMDatasets.source(
    data::RasterData{CHELSA2,T},
    future::Projection{S,M};
    layer=first(SimpleSDMDatasets.layers(data, future)),
    timespan=first(SimpleSDMDatasets.timespans(data, future)),
) where {T<:BioClim,S<:CHELSA2Scenario,M<:CHELSA2Model}
    var_code = (layer isa Integer) ? layer : findfirst(isequal(layer), layers(data))

    year_sep = string(timespan.first.value) * "-" * string(timespan.second.value)

    Mstr, Sstr = map(x -> convert(String, split(string(x), ".")[end]), [M, S])

    model_sep = replace(uppercase(Mstr) * "/" * lowercase(Sstr), "_" => "-")

    root = "https://envicloud.wsl.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/$(year_sep)/$(model_sep)/bio/"

    stem = "CHELSA_bio$(var_code)_$(year_sep)_$(lowercase(replace(Mstr, "_" => "-")))_$(lowercase(Sstr))_V.2.1.tif"

    return (
        url=root * stem,
        filename=lowercase(stem),
        outdir=SimpleSDMDatasets.destination(data, future),
    )
end

function destination(
    ::RasterData{P,D}; kwargs...
) where {P<:RasterProvider,D<:RasterDataset}
    Pstr, Dstr = map(x -> convert(String, split(string(x), ".")[end]), [P, D])
    return joinpath(SimpleSDMDatasets._LAYER_PATH, Pstr, Dstr)
end

function SimpleSDMDatasets.destination(
    ::RasterData{P,D}, ::Projection{S,M}; kwargs...
) where {P<:RasterProvider,D<:RasterDataset,S<:FutureScenario,M<:FutureModel}
    Pstr, Dstr, Mstr, Sstr = map(
        x -> convert(String, split(string(x), ".")[end]), [P, D, M, S]
    )
    return joinpath(
        SimpleSDMDatasets._LAYER_PATH,
        Pstr,
        Dstr,
        replace(Sstr, "_" => "-"),
        replace(Mstr, "_" => "-"),
    )
end


function predict_sdm(climate_layers, model, mat)
    layers_to_matrix!(climate_layers, mat)

    I = eachindex(climate_layers[begin].grid)
    pred = EvoTrees.predict(model, mat')

    distribution = SimpleSDMPredictor(
        zeros(Float32, size(climate_layers[begin])); EXTENT...
    )
    distribution.grid[I] = pred[:, 1]

    uncertainty = SimpleSDMPredictor(zeros(Float32, size(climate_layers[begin])); EXTENT...)
    uncertainty.grid[I] = pred[:, 2]

    return rescale(distribution, (0, 1)), rescale(uncertainty, (0, 1))
end

function compute_fit_stats_and_cutoff(distribution, coords, y)
    cutoff = LinRange(extrema(distribution)..., 500)

    filter!(!isnothing, coords)
    coords = convert(Vector{typeof(coords[begin])}, coords)
    
    obs = y .> 0

    tp = zeros(Float64, length(cutoff))
    fp = zeros(Float64, length(cutoff))
    tn = zeros(Float64, length(cutoff))
    fn = zeros(Float64, length(cutoff))

    for (i, c) in enumerate(cutoff)
        prd = distribution[coords] .>= c
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

function fit_sdm(climate_layers, occurrence_layer, gbrt)
    presences, pseudoabsences = get_pres_and_abs(occurrence_layer)

    X, y, coords = get_features_and_labels(presences, pseudoabsences, climate_layers)
 
    filter!(!isnothing, coords)
    Xtrain, Ytrain, Xtest, Ytest = test_train_split(X, y)

    model = fit_evotree(
        brt_params(gbrt); x_train=Xtrain, y_train=Ytrain, x_eval=Xtest, y_eval=Ytest
    )
    return model, coords, y
end

function write_stats(statsdict, path)
    json_string = JSON.json(statsdict)
    open(path, "w") do f
        JSON.print(f, json_string)
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

function fit_whitening(layers, mat)
    layers_to_matrix!(layers, mat)
    return MultivariateStats.fit(Whitening, mat)
end

function decorrelate_chelsa(layers, w, matrix)
    layers_to_matrix!(layers, matrix)
    decorrelated_matrix = MultivariateStats.transform(w, matrix)

    new_layers = []
    for l in 1:length(layers)
        tmp = convert(Float32, similar(layers[begin]))
        tmp.grid .= nothing
        for (i, idx) in enumerate(eachindex(tmp.grid))
            tmp.grid[idx] = decorrelated_matrix[l, i]
        end
        push!(new_layers, tmp)
    end
    return new_layers
end

function get_pres_and_abs(presences)
    background = pseudoabsencemask(WithinRadius, presences; distance=120.0)
    buffer = pseudoabsencemask(WithinRadius, presences; distance=25.0)
    bgmask = background .& (.!buffer)
    bgpoints = SpeciesDistributionToolkit.sample(bgmask, floor(Int, 0.5sum(presences)))
    replace!(bgpoints, false => nothing)
    return presences, bgpoints
end

function convert_occurrence_to_tif!(species, templatelayer)
    templatelayer.grid .= 0

    occurrence_filename = split(species.name, " ")[1] == "Bombus" ? "bees" : "plants"
    df = CSV.read(
        joinpath(datadir(), "public", "occurrence", "$occurrence_filename.csv"), DataFrame
    )
    thisdf = filter(x -> x.species == species.name, df)
    for r in eachrow(thisdf)
        lat, long = Float32.([r.latitude, r.longitude])
        if check_inbounds(templatelayer, lat, long)
            i = SpeciesDistributionToolkit.SimpleSDMLayers._point_to_cartesian(
                templatelayer, (long, lat)
            )
            templatelayer.grid[i] = 1.0
        end
    end
end
function check_inbounds(tmp, lat, long)
    bb = SimpleSDMLayers.boundingbox(tmp)
    return lat > bb[:bottom] && lat < bb[:top] && long > bb[:left] && long < bb[:right]
end

function SimpleSDMLayers._point_to_cartesian(
    layer::T, c; side=:center
) where {T<:SimpleSDMLayer}
    lonside, latside = true, true
    if side == :bottomleft
        lonside, latside = false, false
    end
    if side == :topright
        lonside, latside = true, true
    end
    lon = SimpleSDMLayers._match_longitude(layer, c[1]; lower=lonside)
    lat = SimpleSDMLayers._match_latitude(layer, c[2]; lower=latside)
    isnothing(lon) && return nothing
    isnothing(lat) && return nothing
    return CartesianIndex(lat, lon)
end

# ==================================================================================
#
#
#   old 
#
#
#
#

#=

function make_sdms()
    run(`mkdir -p $(joinpath("/scratch", "mcatchen", "SDMs"))`)

    bees, plants = get_bee_species(), get_plant_species()
    species = vcat(convert.(Vector{String}, ([bees, plants]))...)

    current_layers = [
        SimpleSDMPredictor(RasterData(CHELSA2, BioClim); layer=l, EXTENT...) for
        l in BIOLAYERS
    ]
    I = common_Is(current_layers)
    mat = zeros(Float32, length(BIOLAYERS), length(I))
    w = fit_whitening(current_layers)
    current_decorrelated_layers = decorrelate_chelsa(current_layers, w, mat)

    models = Dict()

    progbar = ProgressMeter.Progress(length(species))
    @info "Fitting SDMs..."
    for sp in species
        run(`mkdir -p $(joinpath("/scratch", "mcatchen", "SDMs", sp))`)

        occ = convert_occurrence_to_tif(sp, current_decorrelated_layers[begin])
        pres, absen = get_pres_and_abs(occ)
        model, xy, y, xy_pres, presences = fit_sdm(pres, absen, current_decorrelated_layers)

        I = CartesianIndices(size(pres))
        sdm, uncert = predict_sdm(current_decorrelated_layers, model, I)

        models[sp] = model

        run(`mkdir -p $(joinpath("/scratch", "mcatchen", "SDMs", sp, "current"))`)

        _write_geotiff(
            joinpath("/scratch", "mcatchen", "SDMs", sp, "current", "sdm.tif"), sdm
        )
        _write_geotiff(
            joinpath("/scratch", "mcatchen", "SDMs", sp, "current", "uncertainty.tif"),
            uncert,
        )

        dict, τ = compute_fit_stats_and_cutoff(sdm, xy, y)
        write_stats(
            dict, joinpath("/scratch", "mcatchen", "SDMs", sp, "current", "fit.json")
        )

        ProgressMeter.next!(progbar; showvalues=[(Symbol("Species"), sp)])
    end

    @info "end of fit"
    progbar = ProgressMeter.Progress(length(species) * length(SSPs) * (length(years)))
    for y in years
        for ssp in SSPs
            layers = [
                SimpleSDMPredictor(
                    RasterData(CHELSA2, BioClim),
                    Projection(ssp, GFDL_ESM4);
                    timespan=y,
                    layer=l,
                    EXTENT...,
                ) for l in BIOLAYERS
            ]
            theselayers = decorrelate_chelsa(layers, w, mat)
            for sp in species
                @info "projecting $sp in $y $ssp"

                sdm, uncert = predict_sdm(theselayers, models[sp], I)

                ypath = string(y.first.value, "-", y.second.value)
                ssppath = string(ssp)

                run(`mkdir -p $(joinpath("/scratch", "mcatchen", "SDMs", sp, ypath))`)
                run(
                    `mkdir -p $(joinpath("/scratch", "mcatchen", "SDMs", sp, ypath, ssppath))`,
                )

                sdm_path = joinpath(
                    "/scratch", "mcatchen", "SDMs", sp, ypath, ssppath, "sdm.tif"
                )
                uncert_path = joinpath(
                    "/scratch", "mcatchen", "SDMs", sp, ypath, ssppath, "uncertainty.tif"
                )

                _write_geotiff(sdm_path, sdm)
                _write_geotiff(uncert_path, uncert)

                ProgressMeter.next!(
                    progbar;
                    showvalues=[
                        (Symbol("Species"), species),
                        (Symbol("Year"), y),
                        (Symbol("SSP"), ssp),
                    ],
                )
            end
        end
    end
end

function get_bee_species()
    return unique(
        CSV.read(joinpath(datadir(), "occurrence", "bees.csv"), DataFrame).species
    )
end
function get_plant_species()
    return unique(
        CSV.read(joinpath(datadir(), "occurrence", "plants.csv"), DataFrame).species
    )
end

#  ===================================================
#
#  Fit sdm from from QCBONCaseStudy
#
# =============================================

function fit_sdm(presences, absences, climate_layers)
    presences = mask(presences, climate_layers[begin])
    absences = mask(absences, climate_layers[begin])

    xy_presence = keys(replace(presences, false => nothing))
    xy_absence = keys(replace(absences, false => nothing))
    xy = vcat(xy_presence, xy_absence)

    X = hcat([layer[xy] for layer in climate_layers]...)
    y = vcat(fill(1.0, length(xy_presence)), fill(0.0, length(xy_absence)))

    train_size = floor(Int, 0.7 * length(y))
    train_idx = StatsBase.sample(1:length(y), train_size; replace=false)
    test_idx = setdiff(1:length(y), train_idx)

    Xtrain, Xtest = X[train_idx, :], X[test_idx, :]
    Ytrain, Ytest = y[train_idx], y[test_idx]

    GAUSS_TREE_PARAMS = EvoTreeGaussian(;
        loss=:gaussian,
        metric=:gaussian,
        nrounds=100,
        nbins=100,
        λ=0.0,
        γ=0.0,
        η=0.1,
        max_depth=7,
        min_weight=1.0,
        rowsample=0.5,
        colsample=1.0,
    )
    model = fit_evotree(
        GAUSS_TREE_PARAMS; x_train=Xtrain, y_train=Ytrain, x_eval=Xtest, y_eval=Ytest
    )
    return model, xy, y, xy_presence, presences
end

function predict_sdm(climate_layers, model, I)
    all_values = zeros(Float32, length(I), length(climate_layers))

    for (i, idx) in enumerate(I)
        for l in eachindex(climate_layers)
            all_values[i, l] = climate_layers[l].grid[idx]
        end
    end

    pred = EvoTrees.predict(model, all_values)
    distribution = SimpleSDMPredictor(
        zeros(Float32, size(climate_layers[begin]));
        SpeciesDistributionToolkit.boundingbox(climate_layers[begin])...,
    )
    distribution.grid[I] = pred[:, 1]
    distribution

    uncertainty = SimpleSDMPredictor(
        zeros(Float32, size(climate_layers[begin]));
        SpeciesDistributionToolkit.boundingbox(climate_layers[begin])...,
    )
    uncertainty.grid[I] = pred[:, 2]
    uncertainty

    return rescale(distribution, (0, 1)), rescale(uncertainty, (0, 1))
end

function get_pres_and_abs(presences)
    #absences = rand(SurfaceRangeEnvelope, presences)
    background = pseudoabsencemask(WithinRadius, presences; distance=120.0)
    buffer = pseudoabsencemask(WithinRadius, presences; distance=25.0)
    bgmask = background .& (.!buffer)
    bgpoints = SpeciesDistributionToolkit.sample(bgmask, floor(Int, 0.5sum(presences)))
    replace!(bgpoints, false => nothing)
    return presences, bgpoints
end

function convert_occurrence_to_tif(species, templatelayer)
    filename = split(species, " ")[1] == "Bombus" ? "bees" : "plants"
    df = CSV.read(joinpath(datadir(), "occurrence", "$filename.csv"), DataFrame)
    thisdf = filter(x -> x.species == species, df)
    tmp = similar(templatelayer)
    tmp.grid .= 0
    for r in eachrow(thisdf)
        lat, long = Float32.([r.latitude, r.longitude])
        if check_inbounds(tmp, lat, long)
            i = SimpleSDMLayers._point_to_cartesian(tmp, Point(long, lat))
            tmp.grid[i] = 1.0
        end
    end
    return convert(Bool, tmp)
end

function check_inbounds(tmp, lat, long)
    bb = SimpleSDMLayers.boundingbox(tmp)
    return lat > bb[:bottom] && lat < bb[:top] && long > bb[:left] && long < bb[:right]
end

function decorrelate_chelsa(layers, w, matrix)
    get_matrix_form!(layers, Is, matrix)
    decorrelated_matrix = MultivariateStats.transform(w, matrix)

    new_layers = []
    for l in 1:length(layers)
        tmp = convert(Float32, similar(layers[begin]))
        tmp.grid .= nothing
        tmp.grid[Is] .= decorrelated_matrix[l, :]
        tmp
        push!(new_layers, tmp)
    end
    return new_layers
end

function common_Is(layers)
    Is = []
    for l in layers
        push!(Is, findall(x -> !isnothing(x) && !isnan(x), l.grid))
    end
    return Is = unique(intersect(unique(Is)...))
end

function get_matrix_form!(layers, I, matrix)
    for l in eachindex(layers)
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

function compute_fit_stats_and_cutoff(distribution, coords, y)
    cutoff = LinRange(extrema(distribution)..., 500)

    obs = y .> 0

    tp = zeros(Float64, length(cutoff))
    fp = zeros(Float64, length(cutoff))
    tn = zeros(Float64, length(cutoff))
    fn = zeros(Float64, length(cutoff))

    for (i, c) in enumerate(cutoff)
        prd = distribution[coords] .>= c
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

    return Dict(:rocauc => ROCAUC, :threshold => τ, :J => J[last(findmax(J))]), τ
end

#  ===================================================
#
#  Plot
#
# =============================================
# Acquire data

function makeplt(layer)
    count = GeoJSON.read(read("./data/counties.json", String))
    state = GeoJSON.read(read("./data/states.json", String))
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1];
        xticklabelsvisible=false,
        xticklabelpad=-20.0,
        yticks=34:2:42,
        xticks=-110:2:-102,
        dest="+proj=ortho +lon_0=-105 +lat_0=30",
        lonlims=long,
        latlims=lat,
    )

    l = Matrix{Float32}(layer.grid)

    lats = collect(LinRange(lat[1], lat[2], size(l, 1)))
    longs = collect(LinRange(long[1], long[2], size(l, 2)))

    hm = heatmap!(ga, longs, lats, l')

    poly!(
        ga,
        count;
        strokecolor=:lightgrey,
        linestyle=:dash,
        strokewidth=1.0,
        color=(:grey, 0),
        shading=false,
    )
    poly!(ga, state; strokecolor=:white, strokewidth=2, color=(:blue, 0), shading=false)

    Colorbar(fig[1, end + 1], hm; height=Relative(0.7))
    return current_figure()
end

#futurelayers = [SimpleSDMPredictor(RasterData(CHELSA2, BioClim), Projection(SSP370, GFDL_ESM4); timespan=Year(2071)=>Year(2100), layer=l, bbox...) for l in BIOLAYERS]
#decorrelated_layers = decorrelate_chelsa(layers)
#sp = "Bombus bifarius"
replace!(occ, false => nothing)
heatmap(
    elev;
    colormap = :deep,
    axis = (; aspect = DataAspect()),
    figure = (; resolution = (800, 500)),
)
scatter!(keys(pres); color = :black)
scatter!(keys(absen); color = :red)
current_figure() 
#makeplt(elev)
#f = makeplt(sdm)
# save("$sp.png", f, px_per_unit=3)

# ===========================================================================
#
#
#  copied from src for now
#
#  update: i think this is fixed now that it uses SDT.jl
#
#
# ===========================================================================

function _find_span(n, m, M, pos, side)
    side in [:left, :right, :bottom, :top] ||
        throw(ArgumentError("side must be one of :left, :right, :bottom, top"))

    pos > M && return nothing
    pos < m && return nothing
    stride = (M - m) / n
    centers = (m + 0.5stride):stride:(M - 0.5stride)
    pos_diff = abs.(pos .- centers)
    pos_approx = isapprox.(pos_diff, 0.5stride)
    if any(pos_approx)
        if side in [:left, :bottom]
            span_pos = findlast(pos_approx)
        elseif side in [:right, :top]
            span_pos = findfirst(pos_approx)
        end
    else
        span_pos = last(findmin(abs.(pos .- centers)))
    end
    return (stride, centers[span_pos], span_pos)
end

"""
    geotiff(file, ::Type{LT}; bandnumber::Integer=1, left=nothing, right=nothing, bottom=nothing, top=nothing) where {LT <: SimpleSDMLayer}

The geotiff function reads a geotiff file, and returns it as a matrix of the
correct type. The optional arguments `left`, `right`, `bottom`, and `left` are
defining the bounding box to read from the file. This is particularly useful if
you want to get a small subset from large files.

The first argument is the type of the `SimpleSDMLayer` to be returned.
"""
function _read_geotiff(
    file::AbstractString,
    ::Type{LT};
    bandnumber::Integer = 1,
    left = -180.0,
    right = 180.0,
    bottom = -90.0,
    top = 90.0,
) where {LT <: SimpleSDMLayer}
    try
        ArchGDAL.read(file) do stuff
            wkt = ArchGDAL.importPROJ4(ArchGDAL.getproj(stuff))
            wgs84 = ArchGDAL.importEPSG(4326)
            # The next comparison is complete bullshit but for some reason, ArchGDAL has no
            # mechanism to test the equality of coordinate systems. I sort of understand why,
            # but it's still nonsense. So we are left with checking the string representations.
            if string(wkt) != string(wgs84)
                @warn """The dataset is not in WGS84
                We will convert it to WGS84 using gdal_warp, and write it to a temporary file.
                This is not an apology, this is a warning.
                Proceed with caution.
                """
                newfile = tempname()
                run(
                    `$(GDAL.gdalwarp_path()) $file $newfile -t_srs "+proj=longlat +ellps=WGS84"`,
                )
                file = newfile
            end
        end
    catch err
        @info err
    end

    # This next block is reading the geotiff file, but also making sure that we
    # clip the file correctly to avoid reading more than we need.
    layer = ArchGDAL.read(file) do dataset
        transform = ArchGDAL.getgeotransform(dataset)
        # wkt = ArchGDAL.getproj(dataset)

        # The data we need is pretty much always going to be stored in the first
        # band, but this is not the case for the future WorldClim data.
        band = ArchGDAL.getband(dataset, bandnumber)
        T = ArchGDAL.pixeltype(band)

        # We need to check that the nodatavalue is represented in the correct pixeltype,
        # which is not always the case (cough CHELSA2 cough). If this is the case, trying to
        # convert the nodata value will throw an InexactError, so we can catch it and to
        # something about it.
        nodata = ArchGDAL.getnodatavalue(band)
        nodata = isnothing(nodata) ? typemin(T) : nodata

        # Get the correct latitudes
        minlon = transform[1]
        maxlat = transform[4]
        maxlon = minlon + size(band, 1) * transform[2]
        minlat = maxlat - abs(size(band, 2) * transform[6])

        left = isnothing(left) ? minlon : max(left, minlon)
        right = isnothing(right) ? maxlon : min(right, maxlon)
        bottom = isnothing(bottom) ? minlat : max(bottom, minlat)
        top = isnothing(top) ? maxlat : min(top, maxlat)

        lon_stride, lat_stride = transform[2], transform[6]

        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)

        lon_stride, left_pos, min_width = _find_span(width, minlon, maxlon, left, :left)
        _, right_pos, max_width = _find_span(width, minlon, maxlon, right, :right)
        lat_stride, top_pos, max_height = _find_span(height, minlat, maxlat, top, :top)
        _, bottom_pos, min_height = _find_span(height, minlat, maxlat, bottom, :bottom)

        max_height, min_height = height .- (min_height, max_height) .+ 1

        # We are now ready to initialize a matrix of the correct type.
        buffer =
            Matrix{T}(undef, length(min_width:max_width), length(min_height:max_height))
        ArchGDAL.read!(
            dataset,
            buffer,
            bandnumber,
            min_height:max_height,
            min_width:max_width,
        )
        buffer = convert(Matrix{Union{Nothing, eltype(buffer)}}, rotl90(buffer))
        replace!(buffer, nodata => nothing)
        return LT(
            buffer,
            left_pos - 0.5lon_stride,
            right_pos + 0.5lon_stride,
            bottom_pos - 0.5lat_stride,
            top_pos + 0.5lat_stride,
        )
    end

    return layer
end

"""
    geotiff(file::AbstractString, layer::SimpleSDMPredictor{T}; nodata::T=convert(T, -9999)) where {T <: Number}

Write a single `layer` to a `file`, where the `nodata` field is set to an
arbitrary value.
"""
function _write_geotiff(
    file::AbstractString,
    layer::SimpleSDMPredictor{T};
    nodata::T = convert(T, -9999),
) where {T <: Number}
    array_t = _prepare_layer_for_burnin(layer, nodata)
    width, height = size(array_t)

    # Geotransform
    gt = zeros(Float64, 6)
    gt[1] = layer.left
    gt[2] = 2stride(layer, 1)
    gt[3] = 0.0
    gt[4] = layer.top
    gt[5] = 0.0
    gt[6] = -2stride(layer, 2)

    # Write
    prefix = first(split(last(splitpath(file)), '.'))
    ArchGDAL.create(prefix;
        driver = ArchGDAL.getdriver("MEM"),
        width = width, height = height,
        nbands = 1, dtype = T,
        options = ["COMPRESS=LZW"]) do dataset
        band = ArchGDAL.getband(dataset, 1)

        # Write data to band
        ArchGDAL.write!(band, array_t)

        # Write nodata and projection info
        ArchGDAL.setnodatavalue!(band, nodata)
        ArchGDAL.setgeotransform!(dataset, gt)
        ArchGDAL.setproj!(dataset, "EPSG:4326")

        # Write !
        return ArchGDAL.write(
            dataset,
            file;
            driver = ArchGDAL.getdriver("GTiff"),
            options = ["COMPRESS=LZW"],
        )
    end
    return file
end

function _prepare_layer_for_burnin(
    layer::SimpleSDMPredictor{T},
    nodata::T,
) where {T <: Number}
    array = replace(layer.grid, nothing => convert(Float32, nodata))
    array = convert(Matrix{Float32}, array)
    dtype = eltype(array)
    array_t = reverse(permutedims(array, [2, 1]); dims = 2)
    return array_t
end

"""
    geotiff(file::AbstractString, layers::Vector{SimpleSDMPredictor{T}}; nodata::T=convert(T, -9999)) where {T <: Number}

Stores a series of `layers` in a `file`, where every layer in a band. See
`geotiff` for other options.
"""
function _write_geotiff(
    file::AbstractString,
    layers::Vector{SimpleSDMPredictor{T}};
    nodata::T = convert(T, -9999),
    driver::String = "GTiff",
) where {T <: Number}
    bands = 1:length(layers)
    SimpleSDMLayers._layers_are_compatible(layers)
    width, height = size(_prepare_layer_for_burnin(layers[1], nodata))

    # Geotransform
    gt = zeros(Float64, 6)
    gt[1] = layers[1].left
    gt[2] = 2stride(layers[1], 1)
    gt[3] = 0.0
    gt[4] = layers[1].top
    gt[5] = 0.0
    gt[6] = -2stride(layers[1], 2)

    # Write
    prefix = first(split(last(splitpath(file)), '.'))
    ArchGDAL.create(prefix;
        driver = ArchGDAL.getdriver("MEM"),
        width = width, height = height,
        nbands = length(layers), dtype = T,
        options = ["COMPRESS=LZW"]) do dataset
        for i in 1:length(bands)
            band = ArchGDAL.getband(dataset, i)

            # Write data to band
            ArchGDAL.write!(band, _prepare_layer_for_burnin(layers[i], nodata))

            # Write nodata and projection info
            ArchGDAL.setnodatavalue!(band, nodata)
        end
        ArchGDAL.setgeotransform!(dataset, gt)
        ArchGDAL.setproj!(dataset, "EPSG:4326")

        # Write !
        return ArchGDAL.write(
            dataset,
            file;
            driver = ArchGDAL.getdriver(driver),
            options = ["COMPRESS=LZW"],
        )
    end
    return file
end

"""
    geotiff(file::AbstractString, layer::SimpleSDMResponse{T}; nodata::T=convert(T, -9999)) where {T <: Number}

Write a single `SimpleSDMResponse` layer to a file.
"""
function _write_geotiff(
    file::AbstractString,
    layer::SimpleSDMResponse{T};
    kwargs...,
) where {T <: Number}
    return _write_geotiff(file, convert(SimpleSDMPredictor, layer); kwargs...)
end

"""
    geotiff(file::AbstractString, layers::Vector{SimpleSDMResponse{T}}; nodata::T=convert(T, -9999)) where {T <: Number}

Write a vector of `SimpleSDMResponse` layers to bands in a file.
"""
function _write_geotiff(
    file::AbstractString,
    layers::Vector{SimpleSDMResponse{T}};
    kwargs...,
) where {T <: Number}
    return _write_geotiff(file, convert.(SimpleSDMPredictor, layers); kwargs...)
end

=#
