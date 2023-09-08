
function compute_overlap(timespan::Type{T}, scenario::Type{S}) where {T,S}
    data = load_data()
    bee_species, plants_species = bees(data), plants(data)
    all_species = vcat(bee_species, plants_species);

    binary_prediction, probability_prediction, empirical = get_metaweb()

    M = BipartiteNetwork(bee_species, plants_species, any.(binary_prediction .∪ empirical))

    sdms = Dict([sp=>load_sdm(sp, timespan, scenario) for sp in all_species])

    int_richness_map = similar(sdms[bee_species[first]])
    int_richness_map.grid .= 0

    int_uncertainty_map = similar(sdms[bee_species[first]])
    int_uncertainty_map.grid .= 0

    sdm_uncertainty_map = similar(sdms[bee_species[first]])
    sdm_uncertainty_map.grid .= 0


    for s in all_species
        sdm_uncertainty_map.grid .+= sdms[s].uncertainty.grid
    end

    for b in bee_species, p in plants_species
        bee_sdm, plant_sdm = sdms[b], sdms[p]

        expected_cooc = bee_sdm.probability.grid .* plant_sdm.probability.grid
        if M[b,p] > int_threshold
            int_richness_map.grid .+= expected_cooc
        end
        H = entropy([M[b,p], 1-M[b,p]])
        int_uncertainty_map .+= H .* (bee_sdm.probability.grid .* plant_sdm.probability.grid)
    end

    ProjectedOverlap(int_richness_map, int_uncertainty_map, sdm_uncertainty_map)
end
