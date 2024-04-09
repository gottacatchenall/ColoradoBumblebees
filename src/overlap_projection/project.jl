
# so this script takes each future scenario and computes 
# a single scneario. run as an array batch job on slurm.

# every scenario is measured relative to the amount of expected interactions at
# each site in the baseline scenario. 

# expected interactions are computed as:
# Sum across each pair species that was observed to interact and those that we predict to
# interact, add the product of each pairs SDM probability to get expected
# interaction richness

# The SDM uncertainty is computed as the sum of the SDM uncertainty across all
# species.

# The interaction uncertainty computed as the expected interaction richness for
# _all pairs of species_ weighted by the uncertainty of the interaciton
# prediction, computed as the entropy of [p_interaction, 1 - p_interaction].
# justification is more uncertain if they are likely to be there and we are more
# _unsure_ about interaction. 

# Finally, combined uncertainty which is the sum of both uncertainties, each of
# which is scaled to [0,1] before.

# We want this to save 
#   - associated with each scenario
#       - dataframe of each pair of species and their co-occurrence overlap in this scenario
#       relative to baseline
#       - map of expected interaction richness 
#       - map of interaction uncertainty (coccurrence probability  * entropy of
#         [p_interact, 1-p_interact])
#       - map of sdm uncertainty  (sum of uncertainty for each sdm)

# this should probably be two scripts: the first does the baseline, the second
# loads the baseline and does the rest (this avoids i/o issues in parallel)

function get_metaweb(best_fit_dir)
    Pbin, P = get_predicted_network(best_fit_dir)
    O = get_empirical_network()

    Pbin, P, O
end

function get_predicted_network(batch_fit_dir)
    bf = ColoradoBumblebees.open(batch_fit_dir)
    mean_prediction = mean([f.predictions.prediction for f in bf.fits])

    mean_prediction_df = copy(bf.fits[1].predictions)
    mean_prediction_df.prediction .= mean_prediction
    mean_threshold = mean([f.fit_stats["threshold"] for f in bf.fits])

    data = load_data()
    all_bees, all_plants = bees(data), plants(data)

    binary_metaweb, probabilistic_metaweb = zeros(Bool, length(all_bees), length(all_plants)), zeros(Float32, length(all_bees), length(all_plants))

    for (i,b) in enumerate(all_bees), (j,p) in enumerate(all_plants)
        y = filter(x -> x.bee == string(b) && x.plant == string(p) , mean_prediction_df).prediction[1]

        binary_metaweb[i,j] = y > mean_threshold
        probabilistic_metaweb[i,j] = y
    end
    binary_metaweb, probabilistic_metaweb
end

function get_empirical_network()
    data = load_data()
    all_bees, all_plants = bees(data), plants(data)

    empirical_metaweb = zeros(Bool, length(all_bees), length(all_plants))

    

    for (i,b) in enumerate(all_bees), (j,p) in enumerate(all_plants)
        empirical_metaweb[i,j] = length(ColoradoBumblebees.interactions(data, b,p)) > 0
    end
    empirical_metaweb
end

function compute_overlap(best_dir_path, timespan::Type{T}, scenario::Type{S}) where {T,S}
    data = load_data()
    bee_species, plants_species = bees(data), plants(data)
    bee_species, plants_species = bee_species[sortperm([b.name for b in bee_species])], plants_species[sortperm([p.name for p in plants_species])]

    all_species = vcat(bee_species, plants_species);

    binary_prediction, probability_prediction, empirical = get_metaweb(best_dir_path)

    M = BipartiteNetwork( Matrix{Bool}(any.(binary_prediction .∪ empirical)), [b.name for b in bee_species], [p.name for p in plants_species],)
    P = BipartiteProbabilisticNetwork(probability_prediction, [b.name for b in bee_species], [p.name for p in plants_species],)

    sdms = Dict([sp=>load_sdm(sp, timespan, scenario) for sp in all_species])

    int_richness_map = similar(sdms[bee_species[begin]].probability)
    int_richness_map.grid .= 0

    int_uncertainty_map = similar(sdms[bee_species[begin]].probability)
    int_uncertainty_map.grid .= 0

    sdm_uncertainty_map = similar(sdms[bee_species[begin]].probability)
    sdm_uncertainty_map.grid .= 0


    for s in all_species
        sdm_uncertainty_map.grid .+= sdms[s].uncertainty.grid
    end

    cooc_df = DataFrame(bee=[], plant=[], mean_cooccurrence=[], var_cooccurrence=[])

    for b in bee_species, p in plants_species
        bee_sdm, plant_sdm = sdms[b], sdms[p]

        expected_cooc = bee_sdm.probability.grid .* plant_sdm.probability.grid

        push!(cooc_df.bee, b.name)
        push!(cooc_df.plant, p.name)
        push!(cooc_df.mean_cooccurrence, mean(expected_cooc))
        push!(cooc_df.var_cooccurrence, var(expected_cooc))

        # this uses cutoff. not sure if that is ideal
        if M[b.name,p.name] == 1
            int_richness_map.grid .+= expected_cooc
        end

        H = StatsBase.entropy([P[b.name,p.name], 1-P[b.name,p.name]])
        int_uncertainty_map.grid .+= H * (bee_sdm.probability.grid .* plant_sdm.probability.grid)
    end

    ProjectedOverlap{T,S}(int_richness_map, int_uncertainty_map, sdm_uncertainty_map, cooc_df)
end
