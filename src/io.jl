function load_gmm(input_file)
    results = JSON.parsefile(input_file)
end 

function read_sdm_baseline(sdm_directory)
    baseline_directory = joinpath(sdm_directory, "baseline")
    @info sdm_directory
    return Dict(
        :range => SDMLayer(joinpath(baseline_directory, "range.tif")),
        :uncertainty => SDMLayer(joinpath(baseline_directory, "uncertainty.tif")),
        :presences => SDMLayer(joinpath(baseline_directory, "presences.tif")),
        :absences => SDMLayer(joinpath(baseline_directory, "absences.tif")),
    )
end

function read_sdm_future(futures_directory, ssp)
    ssp_directory = joinpath(futures_directory, ssp)
    year_paths = readdir(ssp_directory)

    ssp_dict = Dict()

    for yr in year_paths
        @info futures_directory, ssp, yr
        ssp_dict[yr] = Dict(
            :range => SDMLayer(joinpath(ssp_directory, yr, "range.tif")),
            :uncertainty => SDMLayer(joinpath(ssp_directory, yr, "uncertainty.tif"))
        )
    end
    return ssp_dict
end

function read_sdm_futures(sdm_directory)
    futures_directory = joinpath(sdm_directory, "future")

    futures = Dict()
    for ssp in readdir(futures_directory)
        futures[ssp] = read_sdm_future(futures_directory, ssp)
    end
    return futures
end

function read_sdms(artifact_directory, species_name)
    sdm_directory = joinpath(artifact_directory, species_name, "SDMs")
    Dict(
        :baseline => read_sdm_baseline(sdm_directory),
        :future => read_sdm_futures(sdm_directory),
        :metrics => open(joinpath(sdm_directory, "metrics.json"), "r") do f
            return JSON.parse(f)
        end,
    )

end

function read_sdms(artifact_directory)
    species_dirs = readdir(artifact_directory)
    return Dict([
        sp=>read_sdms(artifact_directory, sp)
        for sp in species_dirs]
    )
end


