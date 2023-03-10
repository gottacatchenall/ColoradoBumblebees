function compute_overlap(data::BeeData, baseline::Scenario, future::Scenario)
    bee_baseline = SimpleSDMPredictor(get_sdm_path(baseline, a))
    plant_baseline = SimpleSDMPredictor(get_sdm_path(baseline, b))

    bee_future = SimpleSDMPredictor(get_sdm_path(future, a))
    plant_future = SimpleSDMPredictor(get_sdm_path(future, b))

    return sum(bee_future.grid .* plant_future.grid) /
           sum(bee_baseline.grid .* plant_baseline.grid)
end

function compute_overlap(
    data::BeeData, a::Bee, b::Plant, baseline::Scenario, future::Scenario
)
    bee_baseline = SimpleSDMPredictor(get_sdm_path(baseline, a))
    plant_baseline = SimpleSDMPredictor(get_sdm_path(baseline, b))

    bee_future = SimpleSDMPredictor(get_sdm_path(future, a))
    plant_future = SimpleSDMPredictor(get_sdm_path(future, b))

    return sum(bee_future.grid .* plant_future.grid) /
           sum(bee_baseline.grid .* plant_baseline.grid)
end

# overlaps_this_bee = compute_overlap(data, bees(data)[2], plants(data)[99],scenarios()[1], scenarios()[end]) 
