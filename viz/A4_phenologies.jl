using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using SpeciesInteractionNetworks
using ColorBlendModes
using JSON
using Dates

const SDT = SpeciesDistributionToolkit
const AG = SDT.SimpleSDMPolygons.AG


function plot_gmm(fig, slice, results; title="")
    x = results["data"]["x"]
    y = results["data"]["y"]
    best_k = results["best_model"]["n_components"]
    components = results["best_model"]["components"]
    waics = results["model_comparison"]["waics"]
    

    begindoy = (Date(2025, 5, 1) - Date(2025, 1, 1)).value
    enddoy = (Date(2025, 10, 1) - Date(2025, 1, 1)).value

    row, col = slice

    isfirstcol = col == 1
    islastrow = row == 7

    ylabel = isfirstcol ? "Number of Observations" : ""

    # Create visualization
    #fig = Figure(size=(900, 500))
  
    months = 5:9
    doy_range_by_month = [((Date(2025, m, 1) - Date(2025, 1, 1)).value, (lastdayofmonth(Date(2025, m,1)) - Date(2025, 1, 1)).value) for m in months]
    xticks = ([median(doy_range_by_month[i]) for i in eachindex(months)], [monthname(i) for i in months])

    # Plot 2: Best fit with 95% CI
    ax = Axis(
        fig[slice...],
        title=title,
        titlealign=:left,
        titlefont=:bold_italic,
        aspect=1,
        xlabel="",
        ylabel=ylabel,
        xticks = xticks,
        xticksvisible = false,
        xticklabelsvisible = islastrow,
        xticklabelrotation=π/2,
        xgridvisible=false,
        ygridvisible=false,
        #yticksvisible = isfirstcol,
        #yticklabelsvisible = isfirstcol,
    )
    ymax = maximum(y) + (0.05*maximum(y))
    xlims!(begindoy, enddoy)
    ylims!(0, ymax)
    # Generate prediction curves with CI using posterior samples
    x_pred = range(begindoy, enddoy, length=365*2)
    

    cols = [ :grey98, :grey85,]
    for (i,doy_range) in enumerate(doy_range_by_month)
        poly!(ax, Point2f[(doy_range[1], 0), (doy_range[1], ymax), (doy_range[2]+1, ymax), (doy_range[2]+1, 0)], color = (cols[1 + i % 2], 0.2))
    end 

    # Extract posterior samples for all components
    n_posterior_samples = length(components[1]["mu_samples"])
    y_pred_samples = zeros(n_posterior_samples, length(x_pred))
    
    # Generate predictions for each posterior sample
    for s in 1:n_posterior_samples
        for (j, x_val) in enumerate(x_pred)
            # Initialize to zero for this prediction point
            y_pred_samples[s, j] = 0.0
            
            # Sum over all components using sampled parameters
            for comp in components
                μ_s = comp["mu_samples"][s]
                σ_s = comp["sigma_samples"][s]
                A_s = comp["amplitude_samples"][s]
                y_pred_samples[s, j] += A_s * exp(-(x_val - μ_s)^2 / (2 * σ_s^2))
            end
        end
    end
    
    # Compute mean and 95% CI
    y_pred_mean = vec(mean(y_pred_samples, dims=1))
    y_pred_lower = [quantile(y_pred_samples[:, j], 0.025) for j in 1:length(x_pred)]
    y_pred_upper = [quantile(y_pred_samples[:, j], 0.975) for j in 1:length(x_pred)]
    

    

    # Plot 95% CI band
    band!(ax, x_pred, y_pred_lower, y_pred_upper, 
          color=(:dodgerblue, 0.4), 
          label="95% CI"
    )
    
    # Plot data
    scatter!(ax, x, y, label="Data", markersize=7, color=(:black, 0.5))
    
    # Plot individual components using posterior mean by computing mean of component curves
    colors = [:red, :green, :orange, :purple,]
    for (idx, comp) in enumerate(components)
        # Compute component curves for each posterior sample
        y_component_samples = zeros(n_posterior_samples, length(x_pred))
        
        for s in 1:n_posterior_samples
            μ_s = comp["mu_samples"][s]
            σ_s = comp["sigma_samples"][s]
            A_s = comp["amplitude_samples"][s]
            
            for (j, x_val) in enumerate(x_pred)
                y_component_samples[s, j] = A_s * exp(-(x_val - μ_s)^2 / (2 * σ_s^2))
            end
        end
        
        # Take mean across posterior samples
        y_component_mean = vec(mean(y_component_samples, dims=1))
        
        lines!(
            ax, 
            x_pred, 
            y_component_mean, 
            label="Component $(comp["component_id"])",
            linestyle=:dash,
            linewidth=2,
            color=(colors[mod1(idx, length(colors))], 0.6)
        )
    end
    
    # Plot mean prediction
    lines!(
        ax, 
        x_pred, 
        y_pred_mean, 
        label="Posterior Mean",
        linewidth=3,
        color=:dodgerblue
    )
    
    #axislegend(ax, position=:rt)
    return fig
end

spnames = readdir("artifacts")
paths = [joinpath("artifacts", sp, "phenology.json") for sp in spnames]

cidx = CartesianIndices((1:7, 1:7))

mkpath(joinpath("plots", "phenology"))

for fig_id in 1:4
    f = Figure(size=(2000, 2000), fonts = (; regular = "Roboto"))
    for (i, ci) in enumerate(cidx)
        sp_idx = ((fig_id-1) * prod(size(cidx))) + i + 1
        if sp_idx <= length(spnames)
            res = load_gmm(paths[sp_idx])
            plot_gmm(f, (ci[2], ci[1]), res; title=spnames[sp_idx])
        end
    end 
    f
    save(joinpath(".", "plots", "phenology", "phen$fig_id.png"), f)
end

