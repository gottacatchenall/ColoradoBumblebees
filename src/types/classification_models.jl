abstract type ClassificationModel end 

DrWatson.default_prefix(cm::ClassificationModel) = string(typeof(cm))*"{"




Base.@kwdef struct XGBoost <: ClassificationModel 
    num_rounds = 20
    max_depth = 8
    learning_rate = 0.5
    minimum_child_weight = 1.
end
(xgb::XGBoost)() = XGBOOST(; num_round = xgb.num_rounds, max_depth = xgb.max_depth, eta = xgb.learning_rate, min_child_weight=xgb.minimum_child_weight)

Base.@kwdef struct BoostedRegressionTree <: ClassificationModel 
    num_rounds = 10
    learning_rate = 0.1
    max_depth = 5 
    gamma = 1.0 
end 
(brt::BoostedRegressionTree)() = BOOSTED_REGRESSION_TREE(; num_rounds = brt.num_rounds, eta = brt.learning_rate, max_depth=brt.max_depth, gamma = brt.gamma)

Base.@kwdef struct RandomForest <: ClassificationModel 
    max_depth = -1          # -1: any 
    n_subfeatures = -1   # -1: sqrt(num_features), 0: num_features
    n_trees = 50 
    sampling_fraction = 0.7
    feature_importance = :split 
end 
(rf::RandomForest)() = RANDOM_FOREST(max_depth = rf.max_depth, n_subfeatures = rf.n_subfeatures, sampling_fraction = rf.sampling_fraction, feature_importance = rf.feature_importance)


Base.@kwdef struct LogisticRegression <: ClassificationModel
    lambda = 1.0
    penalty = :l2
end 
(lr::LogisticRegression)() = LOGISTIC_MODEL(lambda=lr.lambda, penalty = lr.penalty)

struct Ensemble <: ClassificationModel end 