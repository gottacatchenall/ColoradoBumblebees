using Turing
using MCMCChains
using Distributions
using Statistics
using JSON
using Dates


@model function gaussian_mixture_model(x, y, k)
    μ = Vector{Float64}(undef, k)
    σ = Vector{Float64}(undef, k)
    A = Vector{Float64}(undef, k)

    for i in 1:k
        # Center of component
        μ[i] ~ Normal(mean(x), std(x))
        # Width of component 
        σ[i] ~ truncated(Normal(2, 2), 0.1, Inf)
        # Amplitude
        A[i] ~ truncated(Normal(0, 1), 0, Inf)
    end
    
    # Observation noise
    σ_obs ~ truncated(Normal(0, 0.1), 0, Inf)
    
    # Likelihood
    for j in eachindex(x)
        # Sum of all Gaussian components
        predicted = sum(A[i]* exp(-(x[j] - μ[i])^2 / (2 * σ[i]^2)) for i in 1:k)
        y[j] ~ Normal(predicted, σ_obs)
    end
end

# Function to compute model comparison metrics
function compare_models(x, y, k_max; n_samples=1000)
    results = Dict()
    
    for k in 1:k_max
        println("Fitting model with k=$k Gaussians...")
        
        model = gaussian_mixture_model(x, y, k)
        chain = sample(model, NUTS(), n_samples)
        
        # Compute WAIC (Widely Applicable Information Criterion)
        log_likelihood = pointwise_loglikelihoods(model, chain)
        waic_val = waic(log_likelihood)
        
        # Store results
        results[k] = Dict(
            "chain" => chain,
            "waic" => waic_val,
            "n_params" => 3 * k + 1  # μ, σ, A for each component + σ_obs
        )        
        println("  WAIC: $(waic_val.waic)")
    end
    
    return results
end

function pointwise_loglikelihoods(model, chain)
    # Extract parameters from chain
    n_samples = size(chain, 1)
    n_data = length(model.args.x)
    k = model.args.k
    
    log_liks = zeros(n_samples, n_data)
    
    # Extract all parameter chains
    μ_chains = [chain[Symbol("μ[$i]")] for i in 1:k]
    σ_chains = [chain[Symbol("σ[$i]")] for i in 1:k]
    A_chains = [chain[Symbol("A[$i]")] for i in 1:k]
    σ_obs_chain = chain[:σ_obs]
    
    for s in 1:n_samples
        # Extract parameters for this sample
        μ_vals = [μ_chains[i][s] for i in 1:k]
        σ_vals = [σ_chains[i][s] for i in 1:k]
        A_vals = [A_chains[i][s] for i in 1:k]
        σ_obs = σ_obs_chain[s]
        
        # Compute log-likelihood for each data point
        for j in 1:n_data
            x_j = model.args.x[j]
            y_j = model.args.y[j]
            predicted = sum(A_vals[i] * exp(-(x_j - μ_vals[i])^2 / (2 * σ_vals[i]^2)) for i in 1:k)
            log_liks[s, j] = logpdf(Normal(predicted, σ_obs), y_j)
        end
    end
    return log_liks
end 

function waic(log_liks)
    _, n_data = size(log_liks)

    # log pointwise predictive density
    lppd = sum(log(mean(exp.(log_liks[:, i]))) for i in 1:n_data)
    # Effective number of parameters
    p_waic = sum(var(log_liks[:, i]) for i in 1:n_data)
    
    # wAIC
    waic_val = -2 * (lppd - p_waic)
    return (waic=waic_val, lppd=lppd, p_waic=p_waic)
end

function fit_gmm(x, y, k_max; n_samples=1_500, burn_in = 500)
    println("Fitting Gaussian Mixture Models...")
    println("="^50)

    m, M = extrema(x)
    x2 = (x .- m) ./ (M - m)
    y2 = y ./ maximum(y)

    all_waics = Dict{Int, Float64}()
    all_models = Dict()
    
    for k in 1:k_max
        println("Fitting k=$k Gaussians...")
        
        model = gaussian_mixture_model(x2, y2, k)
        chain = sample(model, NUTS(), n_samples)

        chain = chain[burn_in:end]
        
        # Compute WAIC
        log_likelihood = pointwise_loglikelihoods(model, chain)
        waic_val = waic(log_likelihood)
        
        all_waics[k] = waic_val.waic
        all_models[k] = chain
        
        println("\tWAIC: $(round(waic_val.waic, digits=3))")
    end
    
    # Find best model
    best_k = argmin(all_waics)
    best_chain = all_models[best_k]
    
    println("\n" * "="^50)
    println("Best model: k=$best_k components")
    println("="^50)
    
    # Extract parameters for best model
    components = []
    for i in 1:best_k
        μ_samples = m .+ (M - m) .* vec(best_chain[Symbol("μ[$i]")])
        σ_samples = (M-m) .* vec(best_chain[Symbol("σ[$i]")])
        A_samples = vec(best_chain[Symbol("A[$i]")])
        
        μ_post = mean(μ_samples)
        σ_post = mean(σ_samples)


        A_post = mean(A_samples)
        

        # Sample posterior for CI computation (take every 10th sample to reduce file size)
        sample_indices = 1:10:length(μ_samples)
        
        component_info = Dict(
            "component_id" => i,
            "mu_mean" => μ_post,
            "mu_ci_lower" => quantile(μ_samples, 0.025),
            "mu_ci_upper" => quantile(μ_samples, 0.975),
            "mu_samples" => collect(μ_samples[sample_indices]),
            "sigma_mean" => σ_post,
            "sigma_ci_lower" => quantile(σ_samples, 0.025),
            "sigma_ci_upper" => quantile(σ_samples, 0.975),
            "sigma_samples" => collect(σ_samples[sample_indices]),
            "amplitude_mean" => A_post,
            "amplitude_ci_lower" => quantile(A_samples, 0.025),
            "amplitude_ci_upper" => quantile(A_samples, 0.975),
            "amplitude_samples" => collect(A_samples[sample_indices])
        )
        
        push!(components, component_info)
    end
    
    # Observation noise
    σ_obs_post = mean(best_chain[:σ_obs])
    
    results = Dict(
        "metadata" => Dict(
            "n_data_points" => length(x),
            "k_max_tested" => k_max,
            "n_mcmc_samples" => n_samples,
            "timestamp" => string(now())
        ),
        "model_comparison" => Dict(
            "waics" => Dict(string(k) => all_waics[k] for k in 1:k_max),
            "best_k" => best_k
        ),
        "best_model" => Dict(
            "n_components" => best_k,
            "observation_noise" => σ_obs_post,
            "components" => components
        ),
        "data" => Dict(
            "x" => x,
            "y" => y2
        )
    )
    return results
end 

function write_gmm(results_dict, output_file)
    # Save to JSON file
    open(output_file, "w") do io
        JSON.print(io, results_dict, 4)
    end
    println("\n" * "="^50)
    println("Results saved to: $output_file")
    println("="^50)
end