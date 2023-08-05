Base.@kwdef struct BoostedRegressionSDM
    loss = :gaussian
    metric = :gaussian
    nrounds = 100
    nbins = 100
    λ = 0.0
    γ = 0.0
    η = 0.1
    max_depth = 7
    min_weight = 1.0
    rowsample = 0.5
    colsample = 1.0
end

function brt_params(brt::BoostedRegressionSDM)
    return EvoTreeGaussian(; (v => getfield(brt, v) for v in fieldnames(typeof(brt)))...)
end
