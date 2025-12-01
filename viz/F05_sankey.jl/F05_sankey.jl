using CairoMakie
#using SpeciesDistributionToolkit
using Statistics
#using SpeciesInteractionNetworks
#using ColorBlendModes
#using JSON
#using Dates
using SankeyMakie
#using Random

# Scenarios and time periods
SCENARIOS = ["SSP126", "SSP245", "SSP370"]
TIME_PERIODS = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

SCENARIO_LABEL = Dict(
    "SSP126" => "SSP 1-2.6", 
    "SSP245" => "SSP 2-4.5",
    "SSP370" => "SSP 3-7.0"
)

range_sizes = CSV.read(joinpath("data", "range_sizes.csv"), DataFrame)

overlaps = CSV.read(joinpath("data", "overlap.csv"), DataFrame)



species_names = unique(range_sizes.species)

color_dict = Dict(
    "Absolute Loser" => :coral,
    "Absolute Winner" => :lightgreen,
    "Relative Winner" => :skyblue,
    "Relative Loser" => :purple
)

function get_winners_and_losers(ssp, timespan)
    baseline_ranges = filter(r-> r.timespan == "baseline", range_sizes)
    future_ranges = filter(r-> r.ssp == ssp && r.timespan == timespan, range_sizes)
    future_overlap = filter(r-> r.ssp == ssp && r.year == timespan, overlaps)

    df = DataFrame(species=[], ssp=[], timespan=[], status=[], relative_range=[], relative_overlap=[])

    for sp in species_names
        future_range_size = filter(x->x.species == sp, future_ranges).range_size[begin]
        baseline_range_size = filter(x->x.species == sp, baseline_ranges).range_size[begin]
    
        overlap_col = occursin("Bombus", sp) ? "bee" : "plant"
        this_species_overlap = filter(r-> r[overlap_col] == sp, future_overlap)
       
        relative_range_change = future_range_size / baseline_range_size
        @info sp
        relative_overlap = median(this_species_overlap.relative_overlap)

        status = ""
        if relative_range_change > 1
            if relative_overlap > 1
                status = "Absolute Winner"
            else # range growth, but lost overlap
                status = "Relative Loser"
            end
        else
            if relative_overlap > 1 # range loss, but increased overlap 
                status = "Relative Winner"
            else # range loss and overlap loss
                status = "Absolute Loser"
            end
        end
        push!(df, (sp, ssp, timespan, status, relative_range_change, relative_overlap))
    end

    return df
end

df = vcat([get_winners_and_losers(ssp, timespan) for ssp in SCENARIOS, timespan in TIME_PERIODS]...)


function prepare_sankey_data(df, ssp_value)
    df_ssp = filter(row -> row.ssp == ssp_value, df)
    
    timespans = sort(unique(df_ssp.timespan))
    all_statuses = sort(unique(df.status))
    
    
    node_id = 1
    node_dict = Dict{Tuple{Any,Any}, Int}()
    node_labels = String[]
    node_status_map = Any[]
    node_counts = Int[]
    
    for timespan in timespans
        for status in all_statuses
            count = nrow(filter(row -> row.timespan == timespan && row.status == status, df_ssp))
            node_dict[(timespan, status)] = node_id
            
            # Only show count if non-zero
            if count > 0
                push!(node_labels, string(count))
            else
                push!(node_labels, "")  
            end
            push!(node_status_map, status)
            push!(node_counts, count)
            node_id += 1
        end
    end
    
    connections = []
    for i in 1:(length(timespans)-1)
        t1 = timespans[i]
        t2 = timespans[i+1]
        
        # Get species at each timespan
        species_at_t1 = filter(row -> row.timespan == t1, df_ssp)
        species_at_t2 = filter(row -> row.timespan == t2, df_ssp)
        
        # Map species to statuses
        species_status_t1 = Dict(row.species => row.status for row in eachrow(species_at_t1))
        species_status_t2 = Dict(row.species => row.status for row in eachrow(species_at_t2))
        
        # Count transitions
        transition_counts = Dict()
        for species in keys(species_status_t1)
            if haskey(species_status_t2, species)
                status1 = species_status_t1[species]
                status2 = species_status_t2[species]
                key = (status1, status2)
                transition_counts[key] = get(transition_counts, key, 0) + 1
            end
        end
        
        for ((status1, status2), count) in transition_counts
            if count > 0
                node1 = node_dict[(t1, status1)]
                node2 = node_dict[(t2, status2)]
                push!(connections, (node1, node2, count))
            end
        end
    end
    
    # Add dummy zero-weight connections to ensure all nodes exist at each timestep 
    for i in 1:(length(timespans)-1)
        t1 = timespans[i]
        t2 = timespans[i+1]
        for status in all_statuses
            node1 = node_dict[(t1, status)]
            node2 = node_dict[(t2, status)]
            if !any(c -> c[1] == node1 && c[2] == node2, connections)
                # Add with weight 0 (will be invisible)
                push!(connections, (node1, node2, 0))
            end
        end
    end
    
    return connections, node_labels, timespans, node_status_map, node_counts
end


begin 

fig = Figure(size=(1400, 1100))
ssp_values = sort(unique(df.ssp))
all_statuses = sort(unique(df.status))
n_ssps = length(ssp_values)

g = GridLayout(fig[1,1])


for (idx, ssp_value) in enumerate(ssp_values)
    col = idx
    connections, labels, timespans, node_status_map, node_counts = prepare_sankey_data(
        filter(x->occursin("Bombus", x.species), df), 
        ssp_value
    )

    ax = Axis(g[1, col], 
        title=SCENARIO_LABEL[ssp_value],
        titlesize=24,
        titlealign=:left,
        alignmode=Outside()
    )
    hidedecorations!(ax)
    hidespines!(ax)
    
    nodecolors = [color_dict[status] for status in node_status_map]
    sankey!(
        ax, 
        connections;
        nodelabels=labels,
        nodecolor=nodecolors,
        linkcolor=SankeyMakie.Gradient(0.4),
        compact=true,
        fontsize=18
    )
end

    for (idx, ssp_value) in enumerate(ssp_values)
        col = idx
        connections, labels, timespans, node_status_map, node_counts = prepare_sankey_data(
            filter(x->!occursin("Bombus", x.species), df), 
            ssp_value
        )
        
        ax = Axis(g[2, col], 
            title=SCENARIO_LABEL[ssp_value],
            titlesize=24,
            titlealign=:left,
            xticks=(1:4, TIME_PERIODS),
            xticklabelrotation=5π/16,
            xgridvisible=false,
            ygridvisible=false,
            xticksvisible=false,
            yticksvisible=false,
            yticklabelsvisible=false,
            xticklabelsize=20,
            alignmode=Outside()
        )
        #hidedecorations!(ax)
        hidespines!(ax)
        
        nodecolors = [color_dict[status] for status in node_status_map]
        sankey!(
            ax, 
            connections;
            nodelabels=labels,
            nodecolor=nodecolors,
            linkcolor=SankeyMakie.Gradient(0.4),
            compact=true,
            fontsize=18
        )
    end  

    Legend(
        g[:,end+1],
        [MarkerElement(color=(c, 0.4), marker=:rect, markersize=40) for c in [:coral, :skyblue, :lightgreen, :purple]],
        ["Absolute Loser", "Relative Winner", "Absolute Winner", "Relative Loser"],
        labelsize=25,
        framevisible=false,
    )


    ax = Axis(g[1,0])
    limits!(ax, -5, 0, 0, 10)
    hidedecorations!(ax)
    hidespines!(ax)
    text!(ax, -4,6, text="Bees", fontsize=38, font=:bold)

    ax = Axis(g[2,0])
    limits!(ax, -5, 0, 0, 10)
    hidedecorations!(ax)
    hidespines!(ax)
    text!(ax, -4.25, 6, text="Plants", fontsize=38, font=:bold)
    
    colsize!(g, 0, Relative(0.12))
    rowgap!(g, 1, Relative(0.05))
    colgap!(g, 2, Relative(0.05))
    colgap!(g, 3, Relative(0.05))

    rowsize!(g, 1, Relative(0.45))

    fig 
end

fig

save(joinpath("plots", "sankey.png"), fig)
