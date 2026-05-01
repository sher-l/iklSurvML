#' Calculate AUC scores of Machine Learning Models in all data
#'
#' @param res.by.ML.Dev.Prog.Sig the results of function ML.Dev.Prog.Sig
#' @param train_data the training data using in ML.Dev.Prog.Sig
#' @param inputmatrix.list A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
#' @param mode Choose MF models: all, single, double
#' @param AUC_time  c(1,2,3,4,5,6,7,······), for 1 year, 2 years, 3 years......We recommend using the shortest survival time among all queues.
#' @param single_ml If the mode is set to "single", you must fill in the following models: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso").
#' @param double_ml1 If the mode is set to "double", you need to fill in the modeling methods here: c('RSF', "StepCox","CoxBoost","Lasso").
#' @param double_ml2 If the mode is set to "double", you need to fill in the modeling methods here: c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
#' @param alpha_for_Enet 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
#' @param direction_for_stepcox  c("both", "backward", "forward")
#' @param auc_cal_method KM, NNE
#'
#' @return AUC Calculation Result
#' @export
#'
#' @examples
#' all.auc.1y <- cal_AUC_ml_res(
#'   res.by.ML.Dev.Prog.Sig = res,
#'   train_data = list_train_vali_Data[["TCGA"]],
#'   inputmatrix.list = list_train_vali_Data, mode = "all",
#'   AUC_time = 1,
#'   auc_cal_method = "KM"
#' )
#'
cal_AUC_ml_res <- function(res.by.ML.Dev.Prog.Sig = NULL, # ML.Dev.Prog.Sig, 函数计算结果
                           train_data, # ML.Dev.Prog.Sig 中的训练集
                           inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
                           mode = NULL, # all, single, double
                           AUC_time = 1, # c(1,2,3,4,5,6,7,······),1年， 2年， 3年····。这里建议使用所有队列中最短的生存时间。
                           single_ml = NULL, # 如果 mode 为single 则必须要填 c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           double_ml1 = NULL, #  如果mode 为double 则这里需要填写建模的 方法。c('RSF', "StepCox","CoxBoost","Lasso")
                           double_ml2 = NULL, # c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           alpha_for_Enet = NULL, # 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
                           direction_for_stepcox = NULL, #  c("both", "backward", "forward")
                           auc_cal_method = "KM" # KM, NNE 默认为KM
) {
  message("--- Data preprocessing ---")
  auc_cal_method <- match.arg(auc_cal_method, c("KM", "NNE"))

  data_names <- names(inputmatrix.list)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(inputmatrix.list))
  }
  normalized_inputs <- lapply(data_names, function(nm) {
    normalize_ml_data_columns(inputmatrix.list[[nm]], paste0("Dataset '", nm, "'"))
  })
  names(normalized_inputs) <- data_names

  sig.gene <- normalize_ml_feature_names(res.by.ML.Dev.Prog.Sig$Sig.genes)
  common_feature <- normalize_ml_feature_names(c("ID", "OS.time", "OS", sig.gene))
  train_data_norm <- normalize_ml_data_columns(train_data, "Training data")
  preprocess_recipe <- res.by.ML.Dev.Prog.Sig$Preprocess.recipe
  if (is.null(preprocess_recipe)) {
    preprocessed_train <- preprocess_train_data(train_data_norm, common_feature, return_recipe = TRUE)
    preprocess_recipe <- preprocessed_train$recipe
  }
  prepped_inputs <- preprocess_data_list(
    normalized_inputs,
    common_feature,
    recipe = preprocess_recipe
  )

  less.os.time <- min(vapply(prepped_inputs, function(x) max(x$OS.time), numeric(1)))
  if (less.os.time <= 365 * AUC_time) {
    stop(paste0(
      "AUC_time is outside the available follow-up window. The shortest ",
      "overall survival time in the provided cohorts is ", less.os.time,
      " days; choose AUC_time < ", round(less.os.time / 365, 3), " years."
    ), call. = FALSE)
  }

  risk_scores <- cal_RS_ml_res(
    res.by.ML.Dev.Prog.Sig = res.by.ML.Dev.Prog.Sig,
    train_data = train_data,
    inputmatrix.list = inputmatrix.list,
    mode = mode,
    single_ml = single_ml,
    double_ml1 = double_ml1,
    double_ml2 = double_ml2,
    alpha_for_Enet = alpha_for_Enet,
    direction_for_stepcox = direction_for_stepcox
  )

  lapply(risk_scores, function(rs.table.list) {
    lapply(rs.table.list, function(x) {
      calculate_survival_roc_from_risk(
        x,
        AUC_time = AUC_time,
        auc_cal_method = auc_cal_method,
        risk_direction = "higher_is_worse"
      )
    })
  })
}
