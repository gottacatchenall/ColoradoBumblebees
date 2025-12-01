
# ============================================================================
# POLYGON HELPERS (county and state lines)
# ============================================================================

function get_polygons()
    us_states = getpolygon(PolygonData(GADM, Countries), country="USA", level=1)
    us_counties = getpolygon(PolygonData(GADM, Countries), country="USA", level=2)
    states_to_include = ["Colorado", "Utah", "Wyoming", "Nebraska", "NewMexico", "Arizona", "Oklahoma", "Texas", "Kansas", "Arizona"]

    bbox = (left=-109.7, right=-101.8, bottom=34.5,top=42.5)
    bbox_poly = SDT.SimpleSDMPolygons._get_polygon_from_bbox(bbox)

    county_polys = intersect(vcat(
        [us_counties["Level One" => s] for s in states_to_include]...
    ), bbox_poly)
    state_polys = intersect(FeatureCollection(
        [us_states[s] for s in states_to_include],
    ), bbox_poly)

    return state_polys, county_polys
end 

# ============================================================================
# RICHNESS HELPERS
# ============================================================================

get_bee_richness(sdms) = sum([v[:baseline][:range] for (s,v) in sdms if occursin("Bombus", s)])
get_plant_richness(sdms) = sum([v[:baseline][:range] for (s,v) in sdms if !occursin("Bombus", s)])
get_total_species_richness(sdms) = sum([v[:baseline][:range] for (s,v) in sdms])

get_bee_uncertainty(sdms) = sum([v[:baseline][:uncertainty] for (s,v) in sdms if occursin("Bombus", s)])
get_plant_uncertainty(sdms) = sum([v[:baseline][:uncertainty] for (s,v) in sdms if !occursin("Bombus", s)])
get_total_uncertainty(sdms) = sum([v[:baseline][:uncertainty] for (s,v) in sdms])


# ============================================================================
# BIVARIATE HELPERS
# ============================================================================


function _palette(; low=colorant"#f2f2f2", high=colorant"#120fe3", breaks=3)
    breakpoints = LinRange(0.0, 1.0, breaks)
    return Makie.ColorSchemes.weighted_color_mean.(breakpoints, high, low)
end

function discretize(layer, n::Integer)
    return (x -> round(Int64, x)).(rescale(layer, 1, n))
end


function make_bivariate(P,U; nbreaks=5, xlabel="", ylabel="", high2=colorant"#759d77", high1=colorant"#6676a1")
    colormap1 = _palette(; high=high1, breaks=nbreaks)
    colormap2 = _palette(; high=high2, breaks=nbreaks)
    colormatrix = [ColorBlendModes.blend.(c1, c2; mode=BlendMultiply) for c1 in colormap1, c2 in colormap2]

    m1 = discretize(quantize(P), nbreaks)
    m2 = discretize(quantize(U), nbreaks)

    category = similar(m1)
    for i in eachindex(category)
        category[i] = LinearIndices(colormatrix)[m1[i], m2[i]]
    end

    
    category = similar(m1)
    for i in eachindex(category)
        category[i] = LinearIndices(colormatrix)[m1[i], m2[i]]
    end

    return category, colormatrix
    #=
    f = Figure(size=(1200, 900))
    g = GridLayout(f[1,1])
    ax = Axis(g[1,1]; aspect=DataAspect(), xticklabelsize = 26, yticklabelsize=26) 
    heatmap!(ax, category, colormap=vec(colormatrix))

    ax_inset = Axis(
        g[1, 2],
        width=Relative(0.8),
        height=Relative(0.8),
        aspect=1,
        halign=0.4,
        valign=0.5,
        xticklabelsvisible=false,
        xticksvisible=false,
        yticklabelsvisible=false,
        yticksvisible=false,
        xlabel = xlabel,
        ylabel = ylabel,
        xlabelsize = 26,
        ylabelsize = 26,
    )
    heatmap!(ax_inset, colormatrix)
    #hidedecorations!(ax_inset)

    colsize!(g, 2, Relative(0.3))

    lines!(ax, county_poly, color=:grey20, linewidth=0.3)
    lines!(ax, state_poly, color=:grey20, linewidth=1.25)
    f=#
end


# ============================================================================
# SDM CHANGE HELPERS
# ============================================================================

function get_sdm_timeseries(sdms, species, ssp)
    base = Bool.(sdms[species][:baseline][:range])

    years = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]
    [base, [Bool.(sdms[species][:future][ssp][y][:range]) for y in years]...]
end


function compute_gains_and_losses(sdm_timeseries)
    change = Int.(similar(sdm_timeseries[begin]))

    base, early, mid1, mid2, late  = sdm_timeseries
    
    gained_early = !base .& early
    gained_mid1 = !base .& mid1 .& !gained_early
    gained_mid2 = !base .& mid2 .& !gained_mid1
    gained_late = !base .& late .& !gained_mid2
    persists = base .& late

    # TODO: this is more complicated than this.
    # Things can be gained and then lost.
    # This is the case for Bombus griseocollis, where there is range obtained
    # early that is then later lost. 
    # The fig _as is_ reports the _first flip_. 
    lost_early = base .& .!early
    lost_mid1 = base .& .!mid1 .& .!lost_early
    lost_mid2 = base .& .!mid2 .& .!lost_mid1
    lost_late = base .& .!late .& .!lost_mid2


    change.grid .= 0

    change.grid[findall(gained_early)] .= 1
    change.grid[findall(gained_mid1)] .= 2
    change.grid[findall(gained_mid2)] .= 3
    change.grid[findall(gained_late)] .= 4
    change.grid[findall(persists)] .= 5
    change.grid[findall(lost_early)] .= 6
    change.grid[findall(lost_mid1)] .= 7
    change.grid[findall(lost_mid2)] .= 8
    change.grid[findall(lost_late)] .= 9

    return change
end
