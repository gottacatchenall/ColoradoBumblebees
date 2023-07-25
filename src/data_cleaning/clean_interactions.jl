"""
    clean_interactions(s::Type{T}) where T<:Site
"""
function clean_interactions(s::Type{T}) where {T<:Site}
    raw_data_path = datadir(joinpath("embargo", "interactions", "raw", _filename(s)))
    @info raw_data_path
    rawdf = _misc_cleanup(s, CSV.read(raw_data_path, DataFrame))
    cleandf = allocate_clean_df(rawdf)
    fill_clean_df!(s, rawdf, cleandf)

    # Filter out species not in occurrence data
    occ_df = CSV.read(joinpath("data", "public", "occurrence", "occurrence.csv"), DataFrame)
    occ_species = unique(occ_df.species)

    filter(s-> s.plant ∈ occ_species && s.pollinator ∈ occ_species, cleandf)

    outdir = datadir(joinpath("embargo", "interactions", "clean"))
    run(`mkdir -p $outdir`)
    return CSV.write(joinpath(outdir, _filename(s)), cleandf)
end

"""
    fill_clean_df!(s::Type{T}, rawdf, cleandf) 
"""
function fill_clean_df!(s::Type{T}, rawdf, cleandf) where {T<:Site}
    for (i, row) in enumerate(eachrow(rawdf))
        cleandf[i, :datetime] = _datetime(s, row)
        cleandf[i, :plant] = _plant(s, row)
        cleandf[i, :pollinator] = _bee(s, row)
        cleandf[i, :elevation] = _elevation(s, row)
        cleandf[i, :latitude] = _latitude(s, row)
        cleandf[i, :longitude] = _longitude(s, row)
    end
end

const sitenames = [
    "Avery.Wash", "Copper.Creek", "Frogs.Pond", "Gothic.Spires", "Rosy.Hills", "Stormy.Pass"
]

"""
    _datetime(::Type{T}, row) where T <: Site

"""
function _datetime(::Type{PikesPeak}, row)
    year, month, day_of_month, time = row.year, row["month.x"], row["day.x"], row.time_field
    return DateTime(Date(year, month, day_of_month), time)
end
function _datetime(::Type{ElkMeadows}, row)
    year, month, day_of_month = row.Year, row.Month, row.Day
    return DateTime(Date(year, month, day_of_month))
end
_datetime(::Type{Gothic}, row) = DateTime(Date(row.year, 01, 01) + Day(row.doy), Time(12))

"""
    _plant(::Type{T}, row) where T <: Site
"""
_plant(::Type{PikesPeak}, row) = row["ack.nam"]
_plant(::Type{ElkMeadows}, row) = row["Plant species name"]
_plant(::Type{Gothic}, row) = begin
    splitplant = split(row["plant.species"], ".")
    return string(splitplant[1], " ", splitplant[2])
end

"""
    _bee(::Type{T}, row) where T <: Site
"""
_bee(::Type{PikesPeak}, row) = row["pol_sp"]
_bee(::Type{ElkMeadows}, row) = row["Insect species name"]
_bee(::Type{Gothic}, row) = string("Bombus ", row.species)

"""
    _elevation(::Type{T}, row) where T <: Site
"""
_elevation(::Type{PikesPeak}, row) = row["ele_m"]
_elevation(::Type{ElkMeadows}, row) = row["ele_m"]
_elevation(::Type{Gothic}, row) = begin
    elevs = [3045, 2940, 2855, 3015, 2894, 2990]
    elevs[findfirst(x -> x == row.site, sitenames)]
end

"""
    _latitude(::Type{T}, row) where T <: Site
"""
_latitude(::Type{PikesPeak}, row) = row["lat"]
_latitude(::Type{ElkMeadows}, row) = row["lat"]
function _latitude(::Type{Gothic}, row)
    lats = [38.97501, 38.95645, 38.94414, 38.96498, 38.92809, 38.99161]
    return lats[findfirst(x -> x == row.site, sitenames)]
end

"""
    _longitude(::Type{T}, row) where T <: Site
"""
_longitude(::Type{PikesPeak}, row) = row["lon"]
_longitude(::Type{ElkMeadows}, row) = row["lon"]
function _longitude(::Type{Gothic}, row)
    longs = [-106.9924, -106.98276, -106.98507, -106.98698, -106.96297, -107.01027]
    return longs[findfirst(x -> x == row.site, sitenames)]
end

"""
    allocate_clean_df(rawdata)
"""
function allocate_clean_df(rawdata)
    return DataFrame(;
        plant=["" for i in 1:nrow(rawdata)],
        pollinator=["" for i in 1:nrow(rawdata)],
        datetime=[DateTime(0) for i in 1:nrow(rawdata)],
        elevation=zeros(Float32, nrow(rawdata)),
        latitude=zeros(Float32, nrow(rawdata)),
        longitude=zeros(Float32, nrow(rawdata)),
    )
end

"""
    _filename(::Type{T}) where T<:Site
"""
_filename(::Type{PikesPeak}) = "pikespeak.csv"
_filename(::Type{Gothic}) = "gothic.csv"
_filename(::Type{ElkMeadows}) = "elkmeadows.csv"

"""
    _misc_cleanup(::Type{T}, rawdata) where T<:Site
"""
function _misc_cleanup(::Type{PikesPeak}, rawdata)
    filter!(r -> !ismissing(r["ack.nam"]), rawdata)
    filter!(r -> split(r["ack.nam"])[2] ∉ ["sp", "sp.", "SP", "NA"], rawdata)
    filter!(r -> split(r["pol_sp"])[2] ∉ ["sp", "sp.", "SP", "NA"], rawdata)

    for r in eachrow(rawdata)
        spl = split(r["ack.nam"], " ")[1:2]
        if occursin("/", spl[2])
            first = string(split(spl[2], "/")[1])
            r["ack.nam"] = string(spl[1], " ", first)
        else
            r["ack.nam"] = string(spl[1], " ", spl[2])
        end
    end
    return rawdata
end

function _misc_cleanup(::Type{ElkMeadows}, rawdata)
    filter!(r -> !ismissing(r["Plant species name"]), rawdata)
    filter!(
        r -> split(r["Plant species name"])[2] ∉ ["sp", "sp.", "UNKNOWN SP", "NA"], rawdata
    )
    filter!(
        r -> split(r["Insect species name"])[2] ∉ ["sp", "sp.", "UNKNOWN SP", "NA"], rawdata
    )

    for r in eachrow(rawdata)
        spl = split(r["Plant species name"], " ")[1:2]
        if occursin("/", spl[2])
            first = string(split(spl[2], "/")[1])
            r["Plant species name"] = string(spl[1], " ", first)
        else
            r["Plant species name"] = string(spl[1], " ", spl[2])
        end

        if r["Insect species name"] == "Bombus (fervidus) californicus"
            r["Insect species name"] = "Bombus californicus"
        end
    end
    return rawdata
end

function _misc_cleanup(::Type{Gothic}, rawdata)
    return rawdata
end
