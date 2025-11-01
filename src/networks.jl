using DataFrames
using CSV
using Dates
using SpeciesInteractionNetworks

function clean_gothic()
    function clean_gothic_plant_name(raw_plant_name)
        genus, sp = String.(split(raw_plant_name, ".")[1:2])
        return genus * " " * sp
    end
    
    function clean_gothic_date(row)
        Date(row.year, 01, 01) + Day(row.doy)
    end
 

    df = CSV.read(joinpath("data", "raw_interactions", "gothic.csv"), DataFrame)
    df.plant_name = clean_gothic_plant_name.(df[!, "plant.species"])
    df.bee_name = ["Bombus " * sp for sp in df.species]
    df.date = [clean_gothic_date(r) for r in eachrow(df)]
    df.dayofyear = [(r.date - Date(r.year, 01, 01)).value for r in eachrow(df)]
    select(
        df,
        [
            :plant_name,
            :bee_name,
            :year,
            :dayofyear,
            :date,
        ]
    )     
end

function clean_pikespeak()
    function clean_pikespeak_plant(raw_plant_name)
        genus, sp = String.(split(raw_plant_name, " ")[1:2])
        return genus * " " * sp
    end
    function clean_pikespeak_date(row)
        year, month, day_of_month = row.year, row["month.x"], row["day.x"]
        Date(year, month, day_of_month)
    end

    df = CSV.read(joinpath("data", "raw_interactions", "pikespeak.csv"), DataFrame)
    dropmissing!(df, "ack.nam")
    df.plant_name = map(clean_pikespeak_plant, df[!, "ack.nam"])
    df.bee_name = df.pol_sp
    df.year = [parse(Int, split(n, "/")[end]) for n in df.date]
    df.date = [clean_pikespeak_date(r) for r in eachrow(df)]

    df.dayofyear = [(r.date - Date(r.year, 01, 01)).value for r in eachrow(df)]

    select(
        df,
        [
            :plant_name,
            :bee_name,
            :year,
            :date,
            :dayofyear
        ]
    )     
end

function clean_elkmeadows()
    df = CSV.read(joinpath("data", "raw_interactions", "elkmeadows.csv"), DataFrame)

    # drop rows with species uncertainty 
    filter!(x->x["Insect species name"] != "Bombus (fervidus) californicus", df)

    function clean_elkmeadows_date(row)
        year, month, day_of_month = row.Year, row.Month, row.Day
        Date(year, month, day_of_month)
    end
    df.plant_name = String.(df[!, "Plant species name"])
    df.bee_name = String.(df[!, "Insect species name"])
    df.year = df.Year
    df.date = [clean_elkmeadows_date(r) for r in eachrow(df)]
    df.dayofyear = [(r.date - Date(r.year, 01, 01)).value for r in eachrow(df)]

    return select(
        df,
        [
            :plant_name,
            :bee_name,
            :year,
            :date,
            :dayofyear,
        ]
    )     
end

function get_interaction_df()
    name_subs = [
        "Cirsium scariosum/scopulorum" => "Cirsium scariosum_scopulorum",
        "Mentzelia multiflora/speciosa" => "Mentzelia multiflora_speciosa",
    ]
    

    df = vcat(
        clean_elkmeadows(),
        clean_pikespeak(),
        clean_gothic()
    )
    
    gbif_plant_names = df.plant_name
    for (s,t) in name_subs
        gbif_plant_names = [replace(p, s=>t) for p in gbif_plant_names]
    end
    df.plant_name = gbif_plant_names
    return df
end 

function get_metaweb_df()
    df = get_interaction_df()

    bees = unique(df.bee_name)
    plants = unique(df.plant_name)
    combos = DataFrame(
        bee_name = repeat(bees, inner=length(plants)),
        plant_name = repeat(plants, outer=length(bees))
    )
    combos.interacts = [(row.bee_name, row.plant_name) in zip(df.bee_name, df.plant_name) for row in eachrow(combos)]
    return combos   
end 

function interacts(plant_species, bee_species)
    occursin("Bombus", plant_species) && throw(ArgumentError("First species should be a plant"))
    occursin("Bombus", bee_species) || throw(ArgumentError("Second species should be a bee"))
    interacts(get_metaweb_df(), plant_species, bee_species)
end

function interacts(metaweb_df, plant_species, bee_species)
    occursin("Bombus", plant_species) && throw(ArgumentError("First species should be a plant"))
    occursin("Bombus", bee_species) || throw(ArgumentError("Second species should be a bee"))
    idx = findfirst(r-> r.plant_name == plant_species && r.bee_name == bee_species, eachrow(metaweb_df))
    return metaweb_df.interacts[idx]
end

function get_species_list(data_dir)
    gbif_df = get_gbif_data(data_dir)
    taxa_df = get_taxa_df(data_dir)
    species_names = []
    for key in unique(gbif_df.speciesKey)
        idx = findfirst(r->r.species_key == key, eachrow(taxa_df))
        if !isnothing(idx)
            push!(species_names, taxa_df.species_name[idx])
        end
    end
    return species_names
end 

function get_gbif_data(data_dir)
    CSV.read(joinpath(data_dir, "gbif.csv"), DataFrame)
end

function get_taxa_df(data_dir)
    CSV.read(joinpath(data_dir, "taxa.csv"), DataFrame)
end

function get_metaweb()
    mw_df = get_metaweb_df()
    b,p = unique(mw_df.bee_name), unique(mw_df.plant_name)
    adj_mat = [mw_df.interacts[findfirst(x->x.bee_name == bi && x.plant_name == pi, eachrow(mw_df))] for bi in b, pi in p]
    return SpeciesInteractionNetwork(Bipartite(b,p), Binary(adj_mat)) 
end 
