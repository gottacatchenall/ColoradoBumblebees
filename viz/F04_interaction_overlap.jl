"""
    This script plots the median overlap in species ranges aggregated into box plots for different SSPs and timeperoids (left panel) and shows the SSP370 change in overlap for each bee species (right) for each time point (panels).
"""

using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates
using Random
using FileIO

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG

include(joinpath("..", "src", "io.jl"))
include(joinpath("..", "src", "networks.jl"))
include("shared.jl")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Scenarios and time periods
SCENARIOS = ["SSP126", "SSP245", "SSP370"]
TIME_PERIODS = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

# Color configuration
SCENARIO_COLORS = [
    colorant"#EE977A",
    colorant"#F5745A", 
    colorant"#E44F3E",
]
RANGE_SIZE_COLORMAP = Makie.ColorSchemes.navia
MAX_RANGE_SIZE = 0.35

# Figure dimensions
FIGURE_SIZE = (2000, 2200)

# Plot styling
AXIS_CONFIG = (
    rightspinevisible = false,
    topspinevisible = false,
    xgridvisible=false,
    ygridvisible=false,
    titlealign = :left,
)

JITTER_SEED = 1

TITLE_SIZE = 40 
RAINCLOUD_TITLE_SIZE = 30
BOMBUS_NAME_SIZE = 30
SUMMARY_YLABEL_SIZE = 34

RAINCLOUD_YLABEL_SIZE = 24
RAINCLOUD_YTICKLABEL_SIZE = 22

TICK_LABEL_SIZE = 34
LEGEND_SIZE = 45

BOXPLOT_WIDTH = 0.1
BOXPLOT_OFFSET = 0.7 
MARKER_SIZE = 22

YLABEL = "Range Overlap Relative to Baseline"

RAINCLOUD_CONFIG = (
    plot_boxplots = true,
    clouds = nothing,
    cloud_width = 0.5,
    markersize = 18,
    side_nudge = 0.35,
    jitter_width = 0.15,
    gap = 0.,
)

REFERENCE_LINE_CONFIG = (
    linewidth = 4,
    linestyle = :dash,
    color = :grey30,
)

# ============================================================================
# DATA PROCESSING FUNCTIONS
# ============================================================================

function compute_overlap_dataframe(sdms)
    interaction_df = filter(x->x.interacts, get_metaweb_df())
    taxa_df = get_taxa_df("./data")
    #filter!(r->r.plant_name ∈ taxa_df.species_name, interaction_df)
    #filter!(r->r.bee_name ∈ taxa_df.species_name, interaction_df)
    filter!(r->r.plant_name ∈ keys(sdms), interaction_df)
    filter!(r->r.bee_name ∈ keys(sdms), interaction_df)

    df = DataFrame(bee=[], plant=[], year=[], ssp=[], relative_overlap=[])
    for ssp in SCENARIOS
        for row in eachrow(interaction_df)
            plant, bee = row.plant_name, row.bee_name
            
            plant_baseline = Bool.(sdms[plant][:baseline][:range])
            bee_baseline = Bool.(sdms[bee][:baseline][:range])

            baseline_overlap = length(findall(bee_baseline .& plant_baseline))
            for yr in TIME_PERIODS

                plant_future = Bool.(sdms[plant][:future][ssp][yr][:range])
                bee_future = Bool.(sdms[plant][:future][ssp][yr][:range])

                rel_overlap = length(findall(bee_future .& plant_future)) / baseline_overlap

                rel_overlap = isinf(rel_overlap) ? NaN : rel_overlap

                push!(df.plant, plant)
                push!(df.bee, bee)
                push!(df.year, yr)
                push!(df.ssp, ssp)
                push!(df.relative_overlap, rel_overlap)
            end
        end
    end
    return df 
end

function calculate_range_size(sdms, species_name, scenario_key, time_period=nothing)
    if isnothing(time_period)
        range_data = Bool.(sdms[species_name][scenario_key][:range])
    else
        range_data = Bool.(sdms[species_name][scenario_key][scenario][time_period][:range])
    end
    return length(findall(range_data)) / length(range_data.grid)
end

function get_baseline_range_size(sdms, species_name)
    baseline_range = Bool.(sdms[species_name][:baseline][:range])
    return length(findall(baseline_range)) / length(baseline_range.grid)
end

function get_future_range_size(sdms, species_name, scenario, time_period)
    future_range = Bool.(sdms[species_name][:future][scenario][time_period][:range])
    return length(findall(future_range)) / length(future_range.grid)
end

function compute_range_size_df(sdms)    
    range_sizes = DataFrame(
        species = [],
        timespan = [],
        ssp = [],
        range_size = []
    )

    for species in keys(sdms)
        range_size = get_baseline_range_size(sdms, species)
        push!(range_sizes, (species, "baseline", "baseline", range_size))
        for ssp in SCENARIOS, timespan in TIME_PERIODS
            range_size = get_future_range_size(sdms, species, ssp, timespan)
            push!(range_sizes, (species, timespan, ssp, range_size))
        end
    end
    return range_sizes
end


function compute_species_order_by_median(overlap_df, species_col, scenario, time_period)
    species_list = unique(overlap_df[!, species_col])
    median_values = map(species_list) do species
        subset = filter(r -> r[species_col] == species && r.ssp == scenario && r.year == time_period, overlap_df)
        median(subset.relative_overlap)
    end
    return species_list[sortperm(median_values)]
end

function prepare_raincloud_data(overlap_df, scenario, time_period, species_order, species_col)
    categories = String[]
    values = Float64[]
    plant_list = String[]

    for species in species_order
        subset = filter(r -> r[species_col] == species && r.ssp == scenario && r.year == time_period, overlap_df)
        for row in eachrow(subset)
            push!(categories, species)
            push!(plant_list, row.plant)
            push!(values, row.relative_overlap)
        end
    end 
    return categories, values, plant_list
end

function compute_median_overlaps(overlap_df, species_list, species_col, scenario, time_period)
    return map(species_list) do species
        subset = filter(r -> r[species_col] == species && r.ssp == scenario && r.year == time_period, overlap_df)
        median(subset.relative_overlap)
    end
end


# ============================================================================
# ANNOTATION HELPER
# ============================================================================

function annotate_boxplot_points!(ax, overlap_df, species_list, species_col, annotations_config)
    for ann_config in annotations_config
        subset = filter(
            r -> r[species_col] == ann_config.species && 
                 r.ssp == ann_config.scenario && 
                 r.year == ann_config.time_period,
            overlap_df
        )
        
        # Calculate position of point
        scenario_idx = findfirst(==(ann_config.scenario), SCENARIOS)
        period_idx = findfirst(==(ann_config.time_period), TIME_PERIODS)
        x_pos = BOXPLOT_OFFSET * period_idx + BOXPLOT_WIDTH * (scenario_idx - 1.75)
        y_pos = median(subset.relative_overlap)

        
        # Annotation styling
        ann_color = get(ann_config, :color, "#444")
        label_offset = get(ann_config, :label_offset, (50, -50))
        arc_height = get(ann_config, :archeight, 0.1)
        arrow_headlength = get(ann_config, :arrowsize, 16)
        ann_style = get(ann_config, :style, Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=arrow_headlength, color=ann_color), ))
        ann_path = get(ann_config, :path, Ann.Paths.Arc(arc_height))
        ann_fontsize = get(ann_config, :fontsize, 26)
        ann_linewidth = get(ann_config, :linewidth, 3)
        ann_shrink = get(ann_config, :shrink, (5, 15))


        # Add the annotation
        annotation!(
            ax,
            label_offset[1], label_offset[2],
            x_pos, y_pos,
            shrink = ann_shrink,
            color = ann_color,
            text = ann_config.species,
            style = ann_style,
            #path = Makie.Ann.Paths.Corner(),
            path = ann_path,
            labelspace = :relative_pixel,
            fontsize = ann_fontsize,
            linewidth = ann_linewidth
        )
    end
end

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

function create_raincloud_axis!(grid_position, time_period, bee_order, overlap_df, sdms, scenario, show_xlabels)
    ax = Axis(
        grid_position;
        AXIS_CONFIG...,
        xticklabelsvisible = show_xlabels,
        xticks = (1.25:length(bee_order)+0.25, bee_order),
        xticklabelrotation = π/2,
        title = time_period,
        titlesize = RAINCLOUD_TITLE_SIZE,
        ylabel = YLABEL,
        ylabelsize = RAINCLOUD_YLABEL_SIZE,
        yticklabelsize = RAINCLOUD_YTICKLABEL_SIZE,
        xticklabelsize = BOMBUS_NAME_SIZE,
    )

    categories, values, plant_names = prepare_raincloud_data(overlap_df, scenario, time_period, bee_order, :bee)
    valid_indices = findall(!isnan, values)
    categories, values, plant_names = categories[valid_indices], values[valid_indices], plant_names[valid_indices]
    
    range_sizes = Dict(species => get_future_range_size(sdms, species, scenario, time_period) for species in unique(plant_names))


    point_colors = [(get(RANGE_SIZE_COLORMAP, (range_sizes[p] / MAX_RANGE_SIZE) - 0.2), 0.45) for p in plant_names]

    limits!(ax, 0.75, length(bee_order) + 0.75, 0, 3.5)
    Random.seed!(JITTER_SEED)
    rainclouds!(ax, categories, Float32.(values); RAINCLOUD_CONFIG..., color = point_colors, plot_boxplots=false)
    hlines!(ax, [1.0]; REFERENCE_LINE_CONFIG...)
    
    return ax
end

function create_density_axis!(grid_position, time_period, bee_order, overlap_df, sdms, scenario, show_xlabels)
    offset_factor = 3.
    
    ax = Axis(
        grid_position;
        xgridvisible=false,
        yticksvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        leftspinevisible=false,
        xlabel = time_period == TIME_PERIODS[end] ? "Relative Range Overlap" : "",
        xlabelsize=32,
        xticklabelsize=24,
        xticks = 0:0.5:3.5,
        yticklabelsize = 22,
        yticks = ([i/offset_factor + 0.025 for i in 1:length(bee_order)], bee_order),
        title = time_period,
        titlealign = :left,
        titlesize=30,
        aspect=1.2,
    )
    
    xlims!(ax, 0, 3)

    categories, values, plant_names = prepare_raincloud_data(overlap_df, scenario, time_period, bee_order, :bee)
    
    valid_indices = findall(!isnan, values)
    categories, values, plant_names = categories[valid_indices], values[valid_indices], plant_names[valid_indices]
    
    unique_cats = unique(categories)
    range_sizes = Dict(species => get_future_range_size(sdms, species, scenario, time_period) for species in bee_order)
    point_colors = [(get(Makie.ColorSchemes.navia, range_sizes[cat] / MAX_RANGE_SIZE), 0.65) for cat in unique_cats]


    for i in reverse(1:length(unique_cats))
        density!(
            ax, 
            values[categories .== unique_cats[i]], 
            offset=i/offset_factor, 
            color=(point_colors[i], 0.5), 
            strokecolor=:grey50, 
            strokewidth=0.25
        )
    end
    vlines!(ax, [1], linestyle=:dash, linewidth=2, color=:grey30)

    return ax
end


function create_summary_axis!(grid_position, title, ylabel, ylim_min, ylim_max)
    ax = Axis(
        grid_position;
        AXIS_CONFIG...,
        xticks = (BOXPLOT_OFFSET:0.67:(BOXPLOT_OFFSET + 0.67 * (length(TIME_PERIODS) - 1)), TIME_PERIODS),
        title = title,
        titlesize = TITLE_SIZE,
        ylabel = ylabel,
        ylabelsize = SUMMARY_YLABEL_SIZE,
        yticklabelsize = TICK_LABEL_SIZE,
        xticklabelsize = TICK_LABEL_SIZE,
    )
    ylims!(ax, ylim_min, ylim_max)
    return ax
end

function add_boxplots!(ax, overlap_df, species_list, species_col)
    for (scenario_idx, scenario) in enumerate(SCENARIOS)
        for (period_idx, time_period) in enumerate(TIME_PERIODS)
            median_values = compute_median_overlaps(overlap_df, species_list, species_col, scenario, time_period)
            
            x_positions = fill(BOXPLOT_OFFSET * period_idx + BOXPLOT_WIDTH * (scenario_idx - 1.75), length(median_values))
            
            boxplot!(
                ax, 
                x_positions,
                median_values,
                width = BOXPLOT_WIDTH,
                color = SCENARIO_COLORS[scenario_idx],
                markersize = MARKER_SIZE,
            )
        end
    end
    hlines!(ax, [1.0]; REFERENCE_LINE_CONFIG...)
end

function add_colorbar!(grid_position)
    #ghost_axis = Axis(grid_position[1, 3])
    #hidedecorations!(ghost_axis)
    #hidespines!(ghost_axis)
    
    Colorbar(
        grid_position[1, :],
        height = 30,
        vertical = false,
        labelsize = 32,
        label = "Range Size",
        flipaxis = false,
        colormap = RANGE_SIZE_COLORMAP,
        alignmode = Outside(),
        ticks = 0:0.05:MAX_RANGE_SIZE,
        ticklabelsize = 26,
        colorrange = (0, MAX_RANGE_SIZE)
    )
end

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Load and process data
sdms = read_sdms("./artifacts")
overlap_df = compute_overlap_dataframe(sdms)

CSV.write("data/overlap.csv", overlap_df)



range_size_df = compute_range_size_df(sdms) 
CSV.write("data/range_sizes.csv", range_size_df)


# Determine bee ordering based on median overlap
bee_order = compute_species_order_by_median(overlap_df, :bee, SCENARIOS[2], TIME_PERIODS[1])
plant_species = unique(overlap_df.plant)

d = Dict(map(plant_species) do species
    subset = filter(r -> r[:plant] == species && r.ssp == "SSP370" && r.year == "2081-2100", overlap_df)
    species=>median(subset.relative_overlap)
end)
sort(d, lt=(a,b)->d[a]>d[b])

# ===============================================
# BOX PLOT ANNOTATIONS
# ===============================================
plant_boxplot_annotations = [
    (
        species = "Asclepias asperula",
        scenario = "SSP370",
        time_period = "2041-2060",
        label = "Asclepias asperula",
        label_offset = (-160, 60),
        archeight=0.4,
    ),
     (
        species = "Asclepias tuberosa",
        scenario = "SSP370",
        time_period = "2081-2100",
        label = "Asclepias tuberosa",
        label_offset = (-200, 40),
        archeight=0.2,
    ),
    (
        species = "Taraxacum officinale",
        scenario = "SSP370",
        time_period = "2081-2100",
        label = "Taraxacum officinale",
        label_offset = (-200, 15),
        archeight=0.05,
    ),
    (
        species = "Ipomopsis aggregata",
        scenario = "SSP245",
        time_period = "2041-2060",
        label = "Ipomopsis aggregata",
        label_offset = (-160, 70),
        archeight=0.4,
    )
]


tuberosa_img = load(joinpath("viz", "plant_images", "tuberosa.png"))
asperula_img = load(joinpath("viz", "plant_images", "asperula.png"))

# ============================================================================
# CREATE VISUALIZATION
# ============================================================================

begin

    fig = Figure(size = FIGURE_SIZE)
    main_grid = GridLayout(fig[1, 1])
    raincloud_grid = GridLayout(main_grid[1, 1])
    summary_grid = GridLayout(main_grid[:, 2])

    # Create raincloud plots for each time period
    raincloud_axes = []
    for (idx, time_period) in enumerate(TIME_PERIODS)
        create_density_axis!(
            raincloud_grid[idx, 1],
            time_period,
            bee_order,
            overlap_df,
            sdms,
            SCENARIOS[2],
            idx == length(TIME_PERIODS),
        )
    end

    # Create bee overlap summary boxplot
    bee_axis = create_summary_axis!(
        summary_grid[1, 1],
        "Bee Overlap with Interacting Species",
        YLABEL,
        0.4, 1.6
    )
    add_boxplots!(bee_axis, overlap_df, bee_order, :bee)

    # Create plant overlap summary boxplot
    plant_axis = create_summary_axis!(
        summary_grid[2, 1],
        "Plant Overlap with Interacting Species",
        YLABEL,
        0, 8.
    )
    add_boxplots!(plant_axis, overlap_df, plant_species, :plant)
    annotate_boxplot_points!(plant_axis, overlap_df, plant_species, :plant, plant_boxplot_annotations)


    tuber_ax = Axis(
        summary_grid[2, 1],
        width=Relative(0.08),
        height=Relative(0.08),
        aspect=DataAspect(),
        halign=0.6,
        valign=0.98,
        spinewidth=1
    )
    hidedecorations!(tuber_ax),
    hidespines!(tuber_ax)
    image!(tuber_ax, rotr90(tuberosa_img))

    asperula_ax = Axis(
        summary_grid[2, 1],
        width=Relative(0.08),
        height=Relative(0.08),
        aspect=DataAspect(),
        halign=0.06,
        valign=0.83,
        spinewidth=1,
    )
    hidedecorations!(asperula_ax)
    hidespines!(asperula_ax)
    image!(asperula_ax, rotr90(asperula_img))


    # Add legend
    Legend(
        summary_grid[3, 1],
        [MarkerElement(color = c, marker = :rect, markersize = 45, strokecolor = :grey50, strokewidth = 1) for c in SCENARIO_COLORS],
        ["\t" * ssp for ssp in SCENARIOS],
        labelsize = LEGEND_SIZE,
        orientation = :horizontal,
        framevisible = false,
        backgroundcolor = :grey94,
        padding = 30,
        colgap = 100,
    )


    Label(
        raincloud_grid[0,0],
        "A",
        justification = :left,
        font = :bold,
        fontsize = 60
    )
    Label(
        summary_grid[0,0],
        "B",
        justification = :left,
        font = :bold,
        fontsize = 60
    )

    # Add colorbar
    colorbar_grid = GridLayout(raincloud_grid[5, 1])
    add_colorbar!(colorbar_grid)

    # Adjust layout spacing
    colsize!(main_grid, 2, Relative(0.6))
    colgap!(main_grid, 1, Relative(0.))
    rowgap!(summary_grid, 2, Relative(0.03))
    rowgap!(summary_grid, 1, Relative(0.02))

    rowgap!(raincloud_grid, 1, Relative(0.04))
    rowgap!(raincloud_grid, 5, Relative(0.047))

    fig
    save(joinpath("plots", "overlap_change.png"), fig)
end





