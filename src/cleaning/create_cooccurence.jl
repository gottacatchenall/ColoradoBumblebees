
# create two SimpleSDMs at 1km aon the extent
# loop through all species 
# 1 if they cocc, 0 if they don't 

function SpeciesDistributionToolkit.SimpleSDMLayers._match_longitude(
    layer::T, lon::K; lower::Bool=true
) where {T<:SimpleSDMLayer,K<:AbstractFloat}
    # lon > 0. && layer.left <= lon <= layer.right || return nothing
    # lon < 0. && layer.left >= lon >= layer.right || return nothing

    lon == layer.left && return 1
    lon == layer.right && return size(layer, 2)
    relative = (lon - layer.left) / (layer.right - layer.left)
    fractional = relative * size(layer, 2) + 1
    if lon in (layer.left):(2stride(layer, 1)):(layer.right)
        if abs(fractional - round(fractional)) < stride(layer, 1)
            fractional = round(fractional)
        end
        f = lower ? floor : ceil
        d = lower ? 1 : 0
        return min(f(Int64, fractional - d), size(layer, 2))
    else
        return min(floor(Int, fractional), size(layer, 2))
    end
end

function create_cooccurence_data()
    template = convert(
        Float32, SimpleSDMPredictor(RasterData(CHELSA2, BioClim); layer="BIO1", EXTENT...)
    )
    occurrence_df = load_occurrence_data()

    allspecies = unique(occurrence_df.species)

    Ibees = findall(x -> contains(x, "Bombus"), allspecies)
    Iplants = findall([i ∉ allspecies[Ibees] for i in unique(occurrence_df.species)])

    bees, plants = allspecies[Ibees], allspecies[Iplants]

    occ1, occ2 = similar(template), similar(template)
    occ1.grid .= 0
    occ2.grid .= 0

    nrows = prod(length.([bees, plants]))

    cooc_df = DataFrame(
        :bee => ["" for i in 1:nrows],
        :plant => ["" for i in 1:nrows],
        :cooccurence => [-1 for i in 1:nrows],
    )

    cursor = 1
    @showprogress for plant in plants
        thisplant = filter(x -> x.species == plant, occurrence_df)
        occ2 = _fill_occ_layer!(occ2, thisplant)

        for bee in bees
            thisbee = filter(x -> x.species == bee, occurrence_df)
            occ1 = _fill_occ_layer!(occ1, thisbee)

            cooc = _cooccurence(occ1, occ2)

            cooc_df.bee[cursor] = bee
            cooc_df.plant[cursor] = plant
            cooc_df.cooccurence[cursor] = cooc
            cursor += 1
        end
    end
    return CSV.write(datadir("embargo", "cooccurence", "clean", "cooccurence.csv"), cooc_df)
end

function _fill_occ_layer!(occ, df)
    occ.grid .= 0
    for r in eachrow(df)
        x = SpeciesDistributionToolkit.SimpleSDMLayers._match_longitude(occ, r.longitude)
        y = SpeciesDistributionToolkit.SimpleSDMLayers._match_latitude(occ, r.latitude)
        # Note: the grid is transposed 
        occ.grid[y, x] = 1.0
    end
    return occ
end

function _cooccurence(occ1, occ2)
    Itrue1 = findall(x -> x == true, occ1.grid)
    Itrue2 = findall(x -> x == true, occ2.grid)
    return length(Itrue1 ∩ Itrue2) > 0
end
