
struct ConfusionMatrix
    tp::Any
    tn::Any
    fp::Any
    fn::Any
    function ConfusionMatrix(a, b, c, d)
        s = a + b + c + d
        return new(a / s, b / s, c / s, d / s)
    end
end

ConfusionMatrix(A::Matrix{Float64}) = ConfusionMatrix(A[1, 1], A[2, 2], A[2, 1], A[1, 2])
ConfusionMatrix(A::Matrix) = ConfusionMatrix(A[1, 1], A[2, 2], A[2, 1], A[1, 2])

tpr(M::ConfusionMatrix) = M.tp / (M.tp + M.fn)
tnr(M::ConfusionMatrix) = M.tn / (M.tn + M.fp)
ppv(M::ConfusionMatrix) = M.tp / (M.tp + M.fp)
npv(M::ConfusionMatrix) = M.tn / (M.tn + M.fn)
fnr(M::ConfusionMatrix) = M.fn / (M.fn + M.tp)
fpr(M::ConfusionMatrix) = M.fp / (M.fp + M.tn)
fdir(M::ConfusionMatrix) = M.fp / (M.fp + M.tp)
fomr(M::ConfusionMatrix) = M.fn / (M.fn + M.tn)
plr(M::ConfusionMatrix) = tpr(M) / fpr(M)
nlr(M::ConfusionMatrix) = fnr(M) / tnr(M)
pt(M::ConfusionMatrix) = sqrt(fpr(M)) / (sqrt(tpr(M)) + sqrt(fpr(M)))
csi(M::ConfusionMatrix) = M.tp / (M.tp + M.fn + M.fp)
prevalence(M::ConfusionMatrix) = (M.tp + M.fn) / (M.tp + M.fn + M.tn + M.fp)
accuracy(M::ConfusionMatrix) = (M.tp + M.tn) / (M.tp + M.tn + M.fp + M.fn)
balanced(M::ConfusionMatrix) = (tpr(M) + tnr(M)) * 0.5
f1(M::ConfusionMatrix) = 2 * (ppv(M) * tpr(M)) / (ppv(M) + tpr(M))
fm(M::ConfusionMatrix) = sqrt(ppv(M) * tpr(M))
informedness(M::ConfusionMatrix) = tpr(M) + tnr(M) - 1.0
markedness(M::ConfusionMatrix) = ppv(M) + npv(M) - 1.0
dor(M::ConfusionMatrix) = plr(M) / nlr(M)
function κ(M::ConfusionMatrix)
    return 2.0 * (M.tp * M.tn - M.fn * M.fp) /
           ((M.tp + M.fp) * (M.fp + M.tn) + (M.tp + M.fn) * (M.fn + M.tn))
end
function mcc(M::ConfusionMatrix)
    return (M.tp * M.tn - M.fp * M.fn) /
           sqrt((M.tp + M.fp) * (M.tp + M.fn) * (M.tn + M.fp) * (M.tn + M.fn))
end

function ∫(x::Array{T}, y::Array{T}) where {T<:Number}
    S = zero(Float64)
    for i in 2:length(x)
        S += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) * 0.5
    end
    return .-S
end

function threshold(obs, pred; levels=500)
    thresholds = LinRange(minimum(pred), maximum(pred), levels)
    M = Vector{ConfusionMatrix}(undef, length(thresholds))
    for (i, τ) in enumerate(thresholds)
        binpred = pred .>= τ
        tp = sum(obs .& binpred)
        tn = sum(.!obs .& .!binpred)
        fp = sum(.!obs .& binpred)
        fn = sum(obs .& .!binpred)
        M[i] = ConfusionMatrix(tp, tn, fp, fn)
    end
    ROCAUC = ∫(fpr.(M), tpr.(M))
    AUPRC = ∫(tpr.(M), ppv.(M))
    𝐌 = M[last(findmax(informedness.(M)))]
    τ = thresholds[last(findmax(informedness.(M)))]
    return 𝐌, ROCAUC, AUPRC, τ
end

function computemeasures(obs, pred; kwargs...)
    M, ROCAUC, AUPRC, optimalthres = threshold(obs, pred; kwargs...)

    meas = Dict(
        :tpr => tpr,
        :tnr => tnr,
        :ppv => ppv,
        :npv => npv,
        :fnr => fnr,
        :fpr => fpr,
        :acc => accuracy,
        :bac => balanced,
        :f1 => f1,
        :mcc => mcc,
        :fm => fm,
        :Y => informedness,
        :mkd => markedness,
        :κ => κ,
        :mcc => mcc,
    )

    results = Dict(:rocauc => ROCAUC, :prauc => AUPRC)
    for (k, v) in meas
        merge!(results, Dict(k => v(M)))
    end
    return results
end

function computemeasures_mlj(pred, obs; kwargs...)
    ypredict = [p.prob_given_ref[2] for p in pred]
    yobs = [x == true for x in obs]

    M, ROCAUC, AUPRC, optimalthres = threshold(yobs, ypredict; kwargs...)

    meas = Dict(
        :tpr => tpr,
        :tnr => tnr,
        :ppv => ppv,
        :npv => npv,
        :fnr => fnr,
        :fpr => fpr,
        :acc => accuracy,
        :bac => balanced,
        :f1 => f1,
        :mcc => mcc,
        :fm => fm,
        :Y => informedness,
        :mkd => markedness,
        :κ => κ,
        :mcc => mcc,
    )

    results = Dict(:rocauc => ROCAUC, :prauc => AUPRC)
    for (k, v) in meas
        merge!(results, Dict(k => v(M)))
    end
    return results[:prauc]
end
