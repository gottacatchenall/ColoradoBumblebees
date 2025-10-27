using JSON
using GBIF
using CSV
using DataFrames

"""
# convert to GBIF order (counterclockwise)
c = [GI.getgeom(GI.getgeom(bbox_poly.geometry, 1), i) for i in 1:5]
AG.toWKT(AG.createpolygon(reverse(c)))
cols = [:red, :green, :blue, :orange, :red]
# POLYGON ((-109.7 34.5,-101.8 34.5,-101.8 42.5,-109.7 42.5,-109.7 34.5))
"""


include(joinpath("src", "networks.jl"))

function load_env()
    dotenv = open(".env", "r") do file
        read(file, String)
    end
    envvars = split(dotenv, "\n")
    
    for ev in envvars
        k, v = split(ev, "=")
        ENV[k] = v
    end
end

function build_taxon_query(taxa)
    ts = GBIF.taxon.(taxa)
    taxon_query = []
    for t in ts
        levels = [:kingdom, :phylum, :class, :order, :family, :genus, :species]
        level = levels[findlast(l -> getfield(t, l) !== missing, levels)]
        push!(taxon_query, String(level) * "Key" => getfield(t, level).second)
    end
    return taxon_query
end

function query_gbif(species_ids)
    #ts = GBIF.taxon.(taxa)
    #species_ids = [t.species[2] for t in ts]

    query_template = JSON.parse(String(read(joinpath("data", "gbif_query_template.json"))))
    query_template["creator"] = ENV["GBIF_USERNAME"]
    

    push!(
        query_template["predicate"]["predicates"], Dict([
            "type" => "in",
            "key" => "TAXON_KEY",
            "values" => species_ids
        ])
    )
    usn, pwd = ENV["GBIF_USERNAME"], ENV["GBIF_PASSWORD"]

    query_path = joinpath(tempdir(), "query.json")


    open(query_path, "w") do io
        JSON.print(io, query_template)
    end 

    cmd = `curl --include --user $(usn):$(pwd) --header "Content-Type: application/json" --data @$(query_path) https://api.gbif.org/v1/occurrence/download/request`
    run(cmd)
end

taxa_df = CSV.read("data/taxa.csv", DataFrame)
load_env()
query_gbif(taxa_df.species_key)

