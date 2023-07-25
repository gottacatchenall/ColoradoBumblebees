function crossvalidation(model, X, y, species_pairs; ens_size=256, batch_size=64,  kwargs...)
    Is = shuffle(1:nrow(X))
    cut = Int32(floor(0.8 * nrow(X)))
    Itrain, Itest = Is[1:cut], Is[(cut + 1):end]
    ytest = [x == true for x in y[Itest]]
    total_ypredict = zeros(length(Is))
    train_ypredict = zeros(length(Itest))
    for i in 1:ens_size
        mach = machine(model, X, y)
        theserows = balance_sample(y, Itrain, batch_size, 0.5)
        fit!(mach; rows=theserows, verbosity=0)
        pred = MLJ.predict(mach; rows=Itest)
        train_ypredict .+= [p.prob_given_ref[2] for p in pred]

        pred = MLJ.predict(mach)
        total_ypredict .+= [p.prob_given_ref[2] for p in pred]
    end
    train_ypredict = train_ypredict ./ (ens_size)

    eval_dict, thres = computemeasures(ytest, train_ypredict)
    eval_dict[:threshold] = thres

    return _eval_dicts_to_csv([eval_dict]), total_ypredict ./ (ens_size)
end


function get_prediction_df(probs, species_pairs, thres)
    species_pairs[!, :probability] .= probs

    species_pairs[!, :thresholded_probability] .= map(x -> x < thres ? 0 : 1, probs)
    return species_pairs
end


function _eval_dicts_to_csv(dicts)
    cols = collect(keys(dicts[begin]))
    df = DataFrame([[] for c in cols], cols)
    for d in dicts
        for (k, v) in d
            push!(df[!, k], v)
        end
    end
    return df
end

function _cv_test_train_split(X, y; proportion=0.8)
    cutoff = Int32(floor(nrow(X) * proportion))
    I = collect(1:nrow(X))
    shuffle!(I)
    return I[(cutoff + 1):end], I[begin:cutoff]
end
