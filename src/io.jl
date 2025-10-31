
function load_sdm(dir_path)
    @info dir_path
    return Dict(
        :prediction => SDMLayer(joinpath(dir_path, "prediction.tif")),
        :uncertainty => SDMLayer(joinpath(dir_path, "uncertainty.tif")),
        :presences => SDMLayer(joinpath(dir_path, "presences.tif")),
        :absences => SDMLayer(joinpath(dir_path, "absences.tif")),
        :metrics => open(joinpath(dir_path, "metrics.json"), "r") do f
            return JSON.parse(f)
        end,
    )
end

function load_sdms(artifact_dir)
    Dict([r=>load_sdm(joinpath(artifact_dir, r)) for r in readdir(artifact_dir)])
end 

function load_gmm(input_file)
    results = JSON.parsefile(input_file)
end 

