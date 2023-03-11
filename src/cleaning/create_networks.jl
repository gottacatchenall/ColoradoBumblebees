
network_id = 0
interaction_id = 0
taxon_id = 0
node_id = 0

function create_interaction_data()
    sites = [PikesPeak, ElkMeadows, Gothic]
    dirpath = datadir("embargo", "interactions", "clean")
    dfs = [CSV.read(datadir(dirpath, _filename(s)), DataFrame) for s in sites]
    _cut_species!.(dfs)

    species_nodes = get_species_nodes(dfs)

    interactions = []

    plantobjs, beeobjs = initialize_plant_and_bee_objects(species_nodes)

    net = initialize_new_network("ColoradoBumblebees")
    for (i, df) in enumerate(dfs)
        for r in eachrow(df)
            plantname, beename, datetime = r.plant, r.pollinator, DateTime(r.datetime)
            plantnode = species_nodes[findfirst(
                node -> node.name == plantname, species_nodes
            )]
            beenode = species_nodes[findfirst(node -> node.name == beename, species_nodes)]
            coord = (r.longitude, r.latitude)
            elev = r.elevation
            push!(
                interactions,
                initialize_new_interaction(
                    net,
                    beeobjs[beenode],
                    plantobjs[plantnode],
                    plantnode,
                    beenode,
                    datetime,
                    coord,
                    sites[i],
                    elev,
                ),
            )
        end
    end

    return collect(values(beeobjs)), collect(values(plantobjs)), interactions
end

function initialize_plant_and_bee_objects(species_nodes)
    plants, bees = Dict(), Dict()

    for s in species_nodes
        isbee = split(s.name, " ")[1] == "Bombus"
        if isbee
            merge!(bees, Dict(s => Bee(s.name, s)))
        else
            merge!(plants, Dict(s => Plant(s.name, s)))
        end
    end
    return plants, bees
end

"""
    get_species_nodes(dfs)
"""
function get_species_nodes(dfs)
    allplants = unique(vcat([convert(Vector{String}, d[!, :plant]) for d in dfs]...))
    allpollinators = unique(
        vcat([convert(Vector{String}, d[!, :pollinator]) for d in dfs]...)
    )
    return initialize_node.(vcat(allplants, allpollinators))
end

"""
    _cut_species!(df)
"""
function _cut_species!(df)
    filter!(x -> x.plant != "Oreochrysum parryi", df)
    filter!(x -> x.plant != "Salix sp", df)
    return filter!(x -> x.pollinator != "Bombus suckleyi", df)
end

"""
    _filename(::Type{T}) where T<:Site
"""
_filename(::Type{PikesPeak}) = "pikespeak.csv"
_filename(::Type{Gothic}) = "gothic_bombus.csv"
_filename(::Type{ElkMeadows}) = "elkmeadows.csv"

function initialize_node(species)
    global node_id
    reftax = initialize_reference_taxon(species)
    node = MangalNode(node_id, species, now(), now(), reftax)
    node_id += 1
    return node
end

function initialize_reference_taxon(speciesname)
    global taxon_id
    reftax = MangalReferenceTaxon(
        taxon_id,
        speciesname,
        Missing(),
        Missing(),
        Missing(),
        Missing(),
        Missing(),
        now(),
        now(),
    )
    taxon_id += 1
    return reftax
end

function initialize_new_network(name)
    global network_id
    net = MangalNetwork(
        network_id,
        false,
        name,
        now(),
        Missing(),
        now(),
        now(),
        Missing(),
        "",
        false,
        Missing(),
    )
    global network_id += 1
    return net
end

function initialize_new_interaction(
    net, bee, plant, plantnode, beenode, datetime, coord, site, elev
)
    global interaction_id
    try 
        thisint = MangalInteraction(
            interaction_id,
            net,
            beenode,
            plantnode,
            datetime,
            coord,
            false,
            :pollination,
            Missing(),
            1,
            now(),
            now(),
            Missing(),
        )
        interaction_id += 1
        return Interaction{site}(bee, plant, thisint, elev, datetime)

    catch
        @info "Failed with $bee, $plant"
        thisint = MangalInteraction(
            interaction_id,
            net,
            beenode,
            plantnode,
            datetime,
            coord,
            false,
            :pollination,
            Missing(),
            1,
            now(),
            now(),
            Missing(),
        )
        return Interaction{site}(bee, plant, thisint, elev, datetime)
    end

end
