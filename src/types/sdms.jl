struct YearRange
    startyear::Year
    endyear::Year
end
function yearranges()
    years = Year.([1970, 2010, 2040, 2070, 2100])
    return [YearRange(years[t - 1] + Year(1), years[t]) for t in 2:length(years)]
end

ssps() = [SSP126, SSP370, SSP585]
ssp_path(::Type{SSP126}) = "SSP126"
ssp_path(::Type{SSP370}) = "SSP370"
ssp_path(::Type{SSP585}) = "SSP585"
ssp_path(::Type{Missing}) = ""

struct Scenario{T}
    years::YearRange
    ssp::Type{T}
end
function scenarios()
    scen = []
    yrs = yearranges()
    push!(scen, Scenario(yrs[1], Missing))
    for s in ssps(), y in yrs[2:end]
        push!(scen, Scenario(y, s))
    end
    return scen
end

function get_output_dir()
    return if contains(run(`hostname`), "narval")
        joinpath("/scratch", "mcatchen", "BeeSDMs") # cluster
    else
        joinpath(datadir("public", "SDMs")) # local
    end
end

function get_sdm_dir(s::Scenario, sp; cluster=false)
    return if cluster
        joinpath(
        "/scratch", "mcatchen", "BeeSDMs", sp.name, year_path(s.years), ssp_path(s.ssp))
    else
        datadir("public", "SDMs", sp.name, year_path(s.years), ssp_path(s.ssp))
    end
end
function get_sdm_path(s::Scenario, sp, ; kwargs...)
    return joinpath(get_sdm_dir(s, sp; kwargs...), "sdm.tif")
end
function get_uncertainty_path(s::Scenario, sp; kwargs...)
    return joinpath(get_sdm_dir(s, sp; kwargs...), "uncertainty.tif")
end
function year_path(y::YearRange)
    y.startyear.value == 1971 && return "current"
    return string(y.startyear.value, "-", y.endyear.value)
end
