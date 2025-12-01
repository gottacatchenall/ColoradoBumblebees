using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates
using Random

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG

include(joinpath("..", "src", "io.jl"))
include(joinpath("..", "src", "networks.jl"))
include("shared.jl")


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

MAX_RELATIVE_OVERLAP = 6

species_names = unique(range_sizes.species)

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


#filter!(x->x.relative_overlap < MAX_RELATIVE_OVERLAP, df)

unique_status = unique(df.status)

[findfirst(isequal(s), unique_status) for s in df.status]

metaweb = get_metaweb()

degrees = degree(metaweb)

max_marker_size = 20
min_marker_size = 7

colors = [:coral, :skyblue, :lightgreen, :gold]

color_dict = Dict(
    "Absolute Loser" => :coral,
    "Absolute Winner" => :lightgreen,
    "Relative Winner" => :skyblue,
    "Relative Loser" => :purple
)


ρ_sum = 0
ct = 0
for (i, yr) in enumerate(TIME_PERIODS)
    for (j, ssp) in enumerate(SCENARIOS)
        this_df = filter(x->x.ssp == ssp && x.timespan == yr, df)
        ρ_sum += cor(this_df.relative_overlap, this_df.relative_range)
        ct += 1
    end
end
@info "Mean correlation: $(ρ_sum/ct)"

begin
    f = Figure(size=(800, 700))
    g = GridLayout(f[1,1])
    gleft = GridLayout(g[1,1])


    for (i, yr) in enumerate(TIME_PERIODS)
        for (j, ssp) in enumerate(SCENARIOS)
            this_df = filter(x->x.ssp == ssp && x.timespan == yr, df)

            g_local = GridLayout(gleft[j,i], rowgap=0, colgap=0)

            istoprow = j == 1
            isbottomrow = j == 3
            isfirstcolumn = i == 1

            ax_top = Axis(
                g_local[1,1],
                topspinevisible=false,
                bottomspinevisible=false,
            )
            hidedecorations!(ax_top)
            hidespines!(ax_top)
            density!(ax_top, this_df.relative_overlap, color=:grey95, strokecolor=:grey50, strokewidth=1)
            vlines!(ax_top, [median(this_df.relative_overlap)], color=:grey70, linestyle=:dot)


            ax_right = Axis(
                g_local[2,2],
            )
            density!(ax_right, this_df.relative_range, direction = :y,  color=:grey95, strokecolor=:grey50, strokewidth=1)
            hlines!(ax_right, [median(this_df.relative_range)], color=:grey70, linestyle=:dot)
            
            hidedecorations!(ax_right)
            hidespines!(ax_right)

            main_ax = Axis(
                g_local[2,1],
                xticklabelsvisible = isbottomrow,
                xticksvisible = isbottomrow,
                yticksvisible = isfirstcolumn,
                yticklabelsvisible = isfirstcolumn,
                ylabel = isfirstcolumn ? "Rel. Range Size" : "",
                xlabel = isbottomrow ? "Rel. Overlap" : "",
                ylabelvisible = isfirstcolumn,
                xticklabelsize=10,
                yticklabelsize=10,
                xgridvisible=false,
                ygridvisible=false
            )
            limits!(main_ax, 0, MAX_RELATIVE_OVERLAP, 0, 1.5)

            poly!(main_ax, [(0,0), (1,0), (1,1), (0,1)], color=(:coral, 0.4))
            poly!(main_ax, [(1,0), (1,1), (7,1), (7,0)], color=(:skyblue, 0.4))
            poly!(main_ax, [(1,1), (1,2), (7,2), (7,1)], color=(:lightgreen, 0.4))
            poly!(main_ax, [(1,1), (1,2), (0,2), (0,1)], color=(:purple, 0.3))
        
            vlines!(main_ax, [1], color=:grey50)
            hlines!(main_ax, [1], color=:grey50)

            markersize = [degrees[x] for x in this_df.species]

            markersize = min_marker_size .+ max_marker_size .*(markersize ./ maximum(markersize))

            col_idx = [colors[findfirst(isequal(s), unique_status) ] for s in this_df.status]

            scatter!(main_ax,  this_df.relative_overlap, this_df.relative_range, markersize=markersize, color=(:grey5, 0.2))
       

            linkyaxes!(main_ax, ax_right)
            linkxaxes!(main_ax, ax_top)

            colsize!(g_local, 1, Relative(0.87))
            rowsize!(g_local, 2, Relative(0.87))
            rowgap!(g_local, 1, Fixed(0.001))
            colgap!(g_local, 1, Fixed(0.001))

            xlims!(ax_top, 0, MAX_RELATIVE_OVERLAP)
            ylims!(ax_right, 0, 1.5)
        end 
    end

    for (i,lab) in enumerate(TIME_PERIODS)
        Label(gleft[0,i], text=lab * "       ", justification=:left)
    end 
    for (i,lab) in enumerate(SCENARIOS)
        Label(gleft[i,0], text=lab, justification=:left)
    end 

    Legend(
        g[2,:],
        [PolyElement(color=(c, 0.4)) for c in [:coral, :skyblue, :lightgreen, :purple]],
        ["Absolute Loser", "Relative Winner", "Absolute Winner", "Relative Loser"],
        orientation=:horizontal,
    )

    for i in 1:4
        colsize!(gleft, i, Relative(0.22))
    end 
    for i in 1:3
        rowsize!(gleft, i, Relative(0.3))
    end 

    f
end
save("plots/winners_and_losers_scatter.png", f)

