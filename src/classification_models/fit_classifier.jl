function batch_fit(
    classifier::ClassificationModel, 
    rep::Union{SpeciesRepresentations,Vector{S}}, 
    embedding_df::DataFrame,
    num_replicates::Integer; 
    kwargs...
) where S<:SpeciesRepresentations
    BatchFit([fit_classifier(classifier, rep, embedding_df; kwargs...) for i in 1:num_replicates]) 
end


function fit_classifier(
    classifier::ClassificationModel, 
    rep::Union{SpeciesRepresentations,Vector{S}}, 
    embedding_df::DataFrame; 
) where S<:SpeciesRepresentations
    model = classifier()
    
    y, X, _ = unpack(embedding_df, ==(:interaction), ∉([:bee, :plant]))
    y = coerce(y, Multiclass{2})

    train_idx, test_idx, catvec = _cv_test_train_split(X)
    predict_df, fit_stats = ensemble_of_balanced_classifiers(model, embedding_df, X, y, train_idx, test_idx, catvec)

    ClassificationFit(classifier, rep, predict_df, fit_stats)
end

function ensemble_of_balanced_classifiers(model, embedding_df, X, y, train_idx, test_idx, catvec; num_models::I=256, batch_size::I=64, train_balance::F=0.5, kwawgs...) where {I<:Integer,F<:AbstractFloat}
    y_bool = [i == true for i in y]
    y_test = Bool[x == true for x in y[test_idx]]
    y_predict_total = zeros(Float32, nrow(X))
    y_predict_test = zeros(Float32, length(y_test))

    mach = machine(model, X, y)
    for _ in 1:num_models

        balanced_train_rows = balance_sample(y_bool, train_idx, batch_size, train_balance)
        fit!(mach; rows=balanced_train_rows, force=true, verbosity=0)

        y_predict = Float32[pdf(i, true) for i in MLJ.predict(mach)]
        y_predict_total += y_predict
        y_predict_test += y_predict[test_idx]
    end

    y_predict_test ./= num_models
    y_predict_total ./= num_models

    stats_dict, thres = compute_fit_stats(y_test, y_predict_test)
    stats_dict[:threshold] = thres


    outdf = copy(embedding_df)

    outdf.prediction = y_predict_total
    outdf.type = catvec
    outdf[!, [:bee, :plant, :interaction, :prediction, :type]], stats_dict
end

function _cv_test_train_split(X; train_proportion=0.7)

    Is = Random.shuffle(Random.seed!(Dates.datetime2unix(now())), 1:nrow(X))
    cut = Int32(floor(train_proportion * nrow(X)))
    Itrain, Itest = Is[1:cut], Is[(cut + 1):end]

    cat = Symbol[:x for i in 1:nrow(X)]

    for i in Itrain
        cat[i] = :train
    end
    for i in Itest
        cat[i] = :test
    end

    return Itrain, Itest, cat
end

