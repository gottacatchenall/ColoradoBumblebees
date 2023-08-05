const EXTENT = (bottom=34.0, top=44.0, left=-110.5, right=-103.5)
const TEMPORAL_INPUT_DIM = 147
const BIOLAYERS = ["BIO$i" for i in 1:19]


const RANDOM_FOREST = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
const XGBOOST = @load XGBoostClassifier verbosity = 0
const BOOSTED_REGRESSION_TREE = @load EvoTreeClassifier pkg = EvoTrees verbosity = 0
const LOGISTIC_MODEL = @load LogisticClassifier pkg = MLJLinearModels verbosity = 0
