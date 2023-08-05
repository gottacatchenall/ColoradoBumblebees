const FIT_STATS_DIR = "fit_stats"
const PREDICTION_DF_DIR = "predictions"

artifactdir() = "./artifacts"

_write_json(writepath, d::Dict) = open(writepath, "w") do f
    write(f, JSON.json(d))
end

# JSON.parsefile does an annoying thing where vectors get saved with each
# element in it's own vector, e.g. [1,2,3,4] get saved as [[1],[2],[3],[4]]. 
# This wrapper for `JSON.parsefile` finds any vectors in the JSON and fixes them
function _read_json(datapath)
    function _loop(dict)
        for (k,v) in dict
            if typeof(v) <: Dict
                dict[k] = _loop(v)
            end
            if typeof(v) <: Vector
                dict[k] = [x[1] for x in v]
            end
        end
        dict
    end
    dict = JSON.parsefile(datapath)
    _loop(dict)
end 



function load(path::String)
    split(path, "/")[end-1] == "classification_fits"  && return _load_classification_fit(path)
    split(path, "/")[end-1] == "species_representations" && return _load_species_representation(path)
end