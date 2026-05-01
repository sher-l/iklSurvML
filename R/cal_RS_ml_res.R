#' Calculate risk scores of Machine Learning Models in all data
#'
#' @param res.by.ML.Dev.Prog.Sig  The results of function ML.Dev.Prog.Sig
#' @param train_data  The training data using in ML.Dev.Prog.Sig
#' @param inputmatrix.list A list contain the data frames (colnames:ID,OS.time,OS,other variables), log2(x+1)， OS.time(day), OS(0/1)
#' @param mode Choose MF models: 'all', 'single', 'double'
#' @param single_ml If the mode is set to "single", you must fill in the following models: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso").
#' @param double_ml1  If the mode is set to "double", you need to fill in the modeling methods here: c('RSF', "StepCox","CoxBoost","Lasso").
#' @param double_ml2   If the mode is set to "double", you need to fill in the modeling methods here: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
#' @param alpha_for_Enet  One of the values from 0.1 to 0.9. c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9). There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'Enet' is not included in the algorithms you choose. 
#' @param direction_for_stepcox  The parameter for the StepCox. One  from "both", "backward", "forward". There are some conditions you could not set this parameter. [1] The mode is 'all'. [2] The mode is 'single' or 'double', but the 'StepCox' is not included in the algorithms you choose. 
#'
#' @return A list of the risk score calculated by the predictive model in each data in the data of the input matrix list.
#' @export
#'
cal_RS_ml_res <- function(res.by.ML.Dev.Prog.Sig = NULL, # ML.Dev.Prog.Sig, 函数计算结果
                          train_data, # ML.Dev.Prog.Sig 中的训练集
                          inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)
                          mode = NULL, # all, single, double
                          single_ml = NULL, # 如果 mode 为single 则必须要填 c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                          double_ml1 = NULL, #  如果 mode 为double 则这里需要填写建模的 方法。c('RSF', "StepCox","CoxBoost","Lasso")
                          double_ml2 = NULL, # c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                          alpha_for_Enet = NULL, # 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
                          direction_for_stepcox = NULL #  c("both", "backward", "forward")
) {
  message("--- Data preprocessing ---")
  single_ml <- normalize_survival_ml_name(single_ml)
  double_ml1 <- normalize_survival_ml_name(double_ml1)
  double_ml2 <- normalize_survival_ml_name(double_ml2)

  if (is.null(res.by.ML.Dev.Prog.Sig) || is.null(res.by.ML.Dev.Prog.Sig$ml.res)) {
    stop("res.by.ML.Dev.Prog.Sig must contain a non-empty ml.res list", call. = FALSE)
  }
  if (is.null(mode) || !mode %in% c("all", "single", "double")) {
    stop("mode must be one of: all, single, double", call. = FALSE)
  }
  requested_models <- select_survival_model_names(
    model_names = names(res.by.ML.Dev.Prog.Sig$ml.res),
    mode = mode,
    single_ml = single_ml,
    double_ml1 = double_ml1,
    double_ml2 = double_ml2,
    alpha_for_enet = alpha_for_Enet,
    direction_for_stepcox = direction_for_stepcox
  )

  data_names <- names(inputmatrix.list)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(inputmatrix.list))
  }
  inputmatrix.list <- lapply(data_names, function(nm) {
    normalize_ml_data_columns(inputmatrix.list[[nm]], paste0("Dataset '", nm, "'"))
  })
  names(inputmatrix.list) <- data_names
  train_data <- normalize_ml_data_columns(train_data, "Training data")

  sig.gene <- normalize_ml_feature_names(res.by.ML.Dev.Prog.Sig$Sig.genes)
  common_feature <- normalize_ml_feature_names(c("ID", "OS.time", "OS", sig.gene))

  preprocess_recipe <- res.by.ML.Dev.Prog.Sig$Preprocess.recipe
  if (is.null(preprocess_recipe)) {
    preprocessed_train <- preprocess_train_data(train_data, common_feature, return_recipe = TRUE)
    train_data <- preprocessed_train$data
    preprocess_recipe <- preprocessed_train$recipe
  } else {
    train_data <- apply_survival_preprocess_recipe(
      train_data,
      recipe = preprocess_recipe,
      label = "Training data",
      common_feature = common_feature
    )
  }
  inputmatrix.list <- preprocess_data_list(
    inputmatrix.list,
    common_feature,
    recipe = preprocess_recipe
  )

  if (!identical(colnames(inputmatrix.list[[1]])[1:3], c("ID", "OS.time", "OS"))) {
    stop("inputmatrix.list datasets must have first columns ID, OS.time, OS", call. = FALSE)
  }

  missing_by_dataset <- vapply(inputmatrix.list, function(x) {
    paste(setdiff(sig.gene, colnames(x)), collapse = ", ")
  }, character(1))
  missing_by_dataset <- missing_by_dataset[nzchar(missing_by_dataset)]
  if (length(missing_by_dataset) > 0) {
    stop(paste0(
      "Some model genes are missing from input cohorts: ",
      paste(paste0(names(missing_by_dataset), " missing ", missing_by_dataset), collapse = "; ")
    ), call. = FALSE)
  }

  val_dd_list <- lapply(inputmatrix.list, function(x) x[, c("OS.time", "OS", sig.gene), drop = FALSE])
  fallback_train_data <- train_data[, c("OS.time", "OS", sig.gene), drop = FALSE]

  riskscore <- list()
  for (model_name in requested_models) {
    rs <- calculate_survival_model_risk_scores(
      model_name = model_name,
      fit = res.by.ML.Dev.Prog.Sig$ml.res[[model_name]],
      val_dd_list = val_dd_list,
      fallback_train_data = fallback_train_data,
      model_info = res.by.ML.Dev.Prog.Sig$Model.info[[model_name]]
    )
    riskscore[[model_name]] <- return_id_to_rs(rs, inputmatrix.list)
  }

  riskscore
}
