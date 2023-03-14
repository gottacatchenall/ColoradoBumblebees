using DrWatson
@quickactivate "ColoradoBumblebees"

using DataFrames, CSV
using Random
using Flux
using Flux: DataLoader
using ProgressMeter

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

phen = phenology(data)

label_df = label_dataframe(data)

ts = [Float32.(v) for (k, v) in phen]


features, labels = [], []

for r in eachrow(label_df)
    push!(features, (Float32.(phen[r.bee]), Float32.(phen[r.plant])))
    push!(labels, Float32(r.interaction))
end


cut = Int32(floor(0.8length(labels)))
Ishuffled = shuffle(1:length(labels))
Itrain, Itest = Ishuffled[begin:cut], Ishuffled[(cut + 1):end]

loader = DataLoader((features[Itrain], labels[Itrain]))


# needs to take both time series and LSTM embed them

loss(xs, ys) = begin    
    pred =  [fullmodel(x)[1] for x in xs]
    Flux.mse(pred, ys)
end
function fullmodel(x) 
    x1,x2 = x
    h = vcat(enc(x1), enc(x2))
    dec(h)
end 

enc = Chain(LSTM(TEMPORAL_INPUT_DIM, 64), Dense(64, 32), Dense(32,8))
dec = Chain(Dense(2 * 8, 1, σ)) #, Dense(32, 1, σ))
opt = ADAM(1e-6)
ps = Flux.params(enc,dec)

nepochs = 300
batch_size = 64

progbar = ProgressMeter.Progress(nepochs)
for e in 1:nepochs
    ord = sample(Itrain, batch_size; replace=false)
    data_batch = (features[ord,:], labels[ord,:])

    Flux.train!(loss, ps, [data_batch], opt)

    ProgressMeter.next!(
        progbar; showvalues=[(Symbol("Batch Loss"),     loss(data_batch...)
        )]
    )
end

predictions, obs = [x[1] for x in fullmodel.(features[Itest])], Bool.(labels[Itest])


# And we pick thresholds in the [0,1] range
thresholds = range(0.0, 1.0; length=500)

# All this is going to be the components of the adjacency matrix at a given threshold
tp = zeros(Float64, length(thresholds))
fp = zeros(Float64, length(thresholds))
tn = zeros(Float64, length(thresholds))
fn = zeros(Float64, length(thresholds))

# Main loop to get the four components
for (i, thr) in enumerate(thresholds)
    pred = vec(predictions .>= thr)
    tp[i] = sum(pred .& obs)
    tn[i] = sum(.!(pred) .& (.!obs))
    fp[i] = sum(pred .& (.!obs))
    fn[i] = sum(.!(pred) .& obs)
end

# Total number of cases
n = tp .+ fp .+ tn .+ fn

# Diagnostic measures
tpr = tp ./ (tp .+ fn)
fpr = fp ./ (fp .+ tn)
tnr = tn ./ (tn .+ fp)
fnr = fn ./ (fn .+ tp)
acc = (tp .+ tn) ./ (n)
racc = ((tn .+ fp) .* (tn .+ fn) .+ (fn .+ tp) .* (fp .+ tp)) ./ (n .* n)
bacc = ((tp ./ (tp .+ fn)) .+ (tn ./ (fp .+ tn))) ./ 2.0
J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0
κ = (acc .- racc) ./ (1.0 .- racc)
threat = tp ./ (tp .+ fn .+ fp)
fomrate = fn ./ (fn .+ tn)
fdirate = fp ./ (fp .+ tp)
ppv = tp ./ (tp .+ fp)
npv = tn ./ (tn .+ fn)

# This bit is here to get the AUC
dx = [reverse(fpr)[i] - reverse(fpr)[i - 1] for i in 2:length(fpr)]
dy = [reverse(tpr)[i] + reverse(tpr)[i - 1] for i in 2:length(tpr)]
AUC = sum(dx .* (dy ./ 2.0))

# Final thresholding results - we pick the value maximizing Youden's J
thr_index = last(findmax(J))
thr_final = thresholds[thr_index]

# Save the validation measures to a plot
validation = Dict{String,Float64}()
validation["ROC-AUC"] = AUC
validation["Threat score"] = threat[thr_index]
validation["Youden's J"] = J[thr_index]
validation["True Positive Rate"] = tpr[thr_index]
validation["True Negative Rate"] = tnr[thr_index]
validation["False Positive Rate"] = fpr[thr_index]
validation["False Negative Rate"] = fnr[thr_index]
validation["Kappa"] = κ[thr_index]
validation["Accuracy"] = acc[thr_index]
validation["Accuracy (random)"] = racc[thr_index]
validation["Accuracy (balanced)"] = bacc[thr_index]
validation["False Discovery Rate"] = fdirate[thr_index]
validation["False Omission Rate"] = fomrate[thr_index]
validation["Positive Predictive Value"] = ppv[thr_index]
validation["Negative Predictive Value"] = npv[thr_index]
