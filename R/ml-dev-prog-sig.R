# Main entry point for developing prognosis signature
#
# This file contains the main function ML.Dev.Prog.Sig

#' Constructing the optimal predictive model for the prognosis with machine learning algorithms
#'
#' @description
#' Develops a robust predictive model using machine learning algorithms
#' to predict patient prognosis from transcriptomic data.
#'
#' @param train_data The training data with ID, OS.time, and OS as the first three columns.
#'   Starting from the fourth column are the variables for model construction.
#'   Expression should be scaled with log2(x+1).
#'   OS.time: survival time in days. OS: survival status (1=Dead, 0=Alive).
#' @param list_train_vali_Data A list containing training and validation data.
#'   Validation data has the same format as training data.
#' @param candidate_genes Character vector of variables for model development.
#'   These should be included in colnames of training data.
#' @param unicox.filter.for.candi Logical. Whether to use univariate Cox regression
#'   to filter prognostic variables. Default is TRUE.
#' @param unicox_p_cutoff P-value threshold for univariate Cox regression. Default is 0.05.
#' @param mode Algorithm mode: "all", "single", or "double".
#'   "all": use all ten algorithms and combinations.
#'   "single": use only one algorithm.
#'   "double": use combination of two algorithms.
#' @param single_ml One of: "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox",
#'   "superpc", "GBM", "survivalsvm", "Ridge", "Lasso".
#' @param alpha_for_Enet Alpha parameter for Enet (0.1 to 0.9).
#' @param direction_for_stepcox Direction for StepCox: "both", "backward", or "forward".
#' @param double_ml1 First algorithm for combination.
#'   Must be one of: "RSF", "StepCox", "CoxBoost", "Lasso".
#' @param double_ml2 Second algorithm for combination.
#'   Must be one of: "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox",
#'   "superpc", "GBM", "survivalsvm", "Ridge", "Lasso".
#' @param nodesize Node size parameter for RSF. Default is 5. Try 5-10.
#' @param seed Random seed for reproducibility.
#' @param cores_for_parallel Number of cores for parallel processing. Default is 6.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{Cindex.res}{ - C-index results for each model in each dataset}
#'     \item{ml.res}{ - Trained model objects}
#'     \item{riskscore}{ - Risk scores for each sample}
#'     \item{Sig.genes}{ - Signature genes used in the model}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with all algorithms
#' result <- ML.Dev.Prog.Sig(
#'   train_data = train,
#'   list_train_vali_Data = data_list,
#'   candidate_genes = gene_list,
#'   mode = "all",
#'   seed = 5201314
#' )
#'
#' # Example with single algorithm
#' result <- ML.Dev.Prog.Sig(
#'   train_data = train,
#'   list_train_vali_Data = data_list,
#'   candidate_genes = gene_list,
#'   mode = "single",
#'   single_ml = "RSF",
#'   seed = 5201314
#' )
#'
#' # Example with double algorithm combination
#' result <- ML.Dev.Prog.Sig(
#'   train_data = train,
#'   list_train_vali_Data = data_list,
#'   candidate_genes = gene_list,
#'   mode = "double",
#'   double_ml1 = "RSF",
#'   double_ml2 = "CoxBoost",
#'   seed = 5201314
#' )
#' }
ML.Dev.Prog.Sig <- function(train_data,
                            list_train_vali_Data,
                            candidate_genes = NULL,
                            unicox.filter.for.candi = NULL,
                            unicox_p_cutoff = NULL,
                            mode = NULL,
                            single_ml = NULL,
                            alpha_for_Enet = NULL,
                            direction_for_stepcox = NULL,
                            double_ml1 = NULL,
                            double_ml2 = NULL,
                            nodesize = NULL,
                            seed = NULL,
                            cores_for_parallel = NULL) {

  # ---- Set default parameters ----
  if (is.null(alpha_for_Enet)) {
    alpha_for_Enet <- 0.1
  }

  if (is.null(cores_for_parallel)) {
    cores_for_parallel <- 6
  }

  if (is.null(direction_for_stepcox)) {
    direction_for_stepcox <- "both"
  }

  if (is.null(unicox_p_cutoff)) {
    unicox_p_cutoff <- 0.05
  }

  if (is.null(unicox.filter.for.candi)) {
    unicox.filter.for.candi <- TRUE
  }

  rf_nodesize <- nodesize

  # ---- Data preprocessing ----

  # Replace '-' with '.' in column names
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    return(x)
  })

  candidate_genes <- gsub("-", ".", candidate_genes)
  colnames(train_data) <- gsub("-", ".", colnames(train_data))

  # Find common features across all datasets
  common_feature <- c("ID", "OS.time", "OS", candidate_genes)

  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  message(paste0("--- the number of raw candidate genes is ", length(candidate_genes), " ---"))
  message(paste0("--- the number of common features is ", length(common_feature) - 3, " ---"))

  # ---- Parameter validation ----
  if (!is.na(rf_nodesize) &
    !is.na(seed) &
    mode %in% c("all", "single", "double") &
    identical(c("ID", "OS.time", "OS"), colnames(train_data)[1:3]) &
    length(candidate_genes) > 0 &
    identical(c("ID", "OS.time", "OS"), common_feature[1:3]) &
    length(common_feature) > 3) {

    message("--- Data preprocessing ---")

    # Preprocess data lists
    list_train_vali_Data <- preprocess_data_list(list_train_vali_Data, common_feature)
    train_data <- preprocess_train_data(train_data, common_feature)

    # ---- Univariate Cox filtering ----
    if (unicox.filter.for.candi) {
      cd.gene <- common_feature[-c(1:3)]
      cd.gene1 <- sig_unicox(
        gene_list = cd.gene,
        inputSet = train_data,
        unicox_pcutoff = unicox_p_cutoff
      )

      message(paste0(
        "--- the number of final unicox filtered candidate genes is ",
        length(cd.gene1), " ---"
      ))
      print(cd.gene1)

      cd.gene2 <- c(common_feature[1:3], cd.gene1)
      common_feature <- cd.gene2
      train_data <- as.data.frame(train_data)[, common_feature]
    } else {
      message(paste0(
        "--- the number of final not unicox filtered candidate genes is ",
        length(common_feature) - 3, " ---"
      ))
    }

    # Prepare estimation and validation data
    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
    val_dd_list <- lapply(list_train_vali_Data, function(x) {
      x[, common_feature[-1]]
    })
    pre_var <- common_feature[-c(1:3)]

    # ---- Run ML algorithms based on mode ----
    if (mode == "all") {
      result <- run_all_algorithms(
        est_dd = est_dd,
        train_data = train_data,
        val_dd_list = val_dd_list,
        list_train_vali_Data = list_train_vali_Data,
        pre_var = pre_var,
        rf_nodesize = rf_nodesize,
        seed = seed,
        cores_for_parallel = cores_for_parallel
      )
    } else if (mode == "single") {
      result <- run_single_algorithm(
        est_dd = est_dd,
        train_data = train_data,
        val_dd_list = val_dd_list,
        list_train_vali_Data = list_train_vali_Data,
        pre_var = pre_var,
        single_ml = single_ml,
        rf_nodesize = rf_nodesize,
        seed = seed,
        alpha_for_enet = alpha_for_Enet,
        direction_for_stepcox = direction_for_stepcox,
        cores_for_parallel = cores_for_parallel
      )
    } else if (mode == "double") {
      result <- run_double_algorithm(
        est_dd = est_dd,
        train_data = train_data,
        val_dd_list = val_dd_list,
        list_train_vali_Data = list_train_vali_Data,
        pre_var = pre_var,
        double_ml1 = double_ml1,
        double_ml2 = double_ml2,
        rf_nodesize = rf_nodesize,
        seed = seed,
        alpha_for_enet = alpha_for_Enet,
        direction_for_stepcox = direction_for_stepcox,
        cores_for_parallel = cores_for_parallel
      )
    }

    message("--- The analysis has been completed ---")
    return(result)
  } else {
    stop(paste(
      "Please provide the full parameters, verify that the column names",
      "match the conditions (ID, OS.time, OS for the first through third columns),",
      "and verify there exist common genes in all cohorts"
    ))
  }
}

#' Run all algorithms and combinations
#'
#' @inheritParams ML.Dev.Prog.Sig
#' @return Result list
#' @keywords internal
run_all_algorithms <- function(est_dd,
                               train_data,
                               val_dd_list,
                               list_train_vali_Data,
                               pre_var,
                               rf_nodesize,
                               seed,
                               cores_for_parallel) {
  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # ---- 1. RSF ----
  message("--- 1-1 RSF ---")
  set.seed(seed)
  fit <- train_rsf(est_dd, rf_nodesize, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_rsf(fit, x))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "RSF")
  result <- rbind(result, cc)
  ml.res[["RSF"]] <- fit
  riskscore[["RSF"]] <- rs

  # ---- RSF combinations ----
  rsf_combos <- c("CoxBoost", "Enet", "GBM", "Lasso", "plsRcox", "Ridge", "StepCox", "SuperPC", "survivalsvm")
  for (algo in rsf_combos) {
    message(paste0("--- 1-", match(algo, rsf_combos) + 1, " RSF + ", algo, " ---"))
    combo_result <- run_rsf_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      rf_nodesize = rf_nodesize,
      seed = seed,
      second_algo = algo,
      cores_for_parallel = cores_for_parallel
    )
    if (!is.null(combo_result)) {
      result <- rbind(result, combo_result$result)
      ml.res <- c(ml.res, combo_result$ml.res)
      riskscore <- c(riskscore, combo_result$riskscore)
    }
  }

  # ---- 2. Enet ----
  message("--- 2. Enet ---")
  for (alpha in seq(0.1, 0.9, 0.1)) {
    message(paste0("--- 2. Enet[α=", alpha, "] ---"))
    set.seed(seed)
    fit <- train_enet(est_dd, pre_var, alpha, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, paste0("Enet[α=", alpha, "]"))
    result <- rbind(result, cc)
    ml.res[[paste0("Enet[α=", alpha, "]")]] <- fit
    riskscore[[paste0("Enet[α=", alpha, "]")]] <- rs
  }

  # ---- 3. StepCox ----
  message("--- 3. StepCox ---")
  for (direction in c("both", "backward", "forward")) {
    message(paste0("--- 3. StepCox[", direction, "] ---"))
    fit <- train_stepcox(est_dd, direction)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(fit, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, paste0("StepCox[", direction, "]"))
    result <- rbind(result, cc)
    ml.res[[paste0("StepCox[", direction, "]")]] <- fit
    riskscore[[paste0("StepCox[", direction, "]")]] <- rs
  }

  # StepCox[both] combinations
  message("--- StepCox[both] combinations ---")
  stepcox_combos <- c("CoxBoost", "Enet", "GBM", "Lasso", "plsRcox", "Ridge", "RSF", "SuperPC", "survivalsvm")
  for (algo in stepcox_combos) {
    combo_result <- run_stepcox_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      direction = "both",
      seed = seed,
      second_algo = algo,
      cores_for_parallel = cores_for_parallel
    )
    if (!is.null(combo_result)) {
      result <- rbind(result, combo_result$result)
      ml.res <- c(ml.res, combo_result$ml.res)
      riskscore <- c(riskscore, combo_result$riskscore)
    }
  }

  # ---- 4-10. Other single algorithms ----
  # CoxBoost
  message("--- 4. CoxBoost ---")
  set.seed(seed)
  fit <- train_coxboost(est_dd, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(fit, x))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "CoxBoost")
  result <- rbind(result, cc)
  ml.res[["CoxBoost"]] <- fit
  riskscore[["CoxBoost"]] <- rs

  # plsRcox
  message("--- 5. plsRcox ---")
  fit <- train_plsrcox(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "plsRcox")
  result <- rbind(result, cc)
  ml.res[["plsRcox"]] <- fit
  riskscore[["plsRcox"]] <- rs

  # SuperPC
  message("--- 6. SuperPC ---")
  superpc_result <- train_superpc(est_dd, seed)
  fit <- superpc_result$fit
  cv_fit <- superpc_result$cv_fit
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = predict_superpc(fit, cv_fit, est_dd, x))
  })
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "SuperPC")
  result <- rbind(result, cc)
  ml.res[["SuperPC"]] <- list(fit, cv_fit)
  riskscore[["SuperPC"]] <- rs

  # GBM
  message("--- 7. GBM ---")
  gbm_result <- train_gbm(est_dd, seed, cores_for_parallel)
  fit <- gbm_result$fit
  best <- gbm_result$best
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(fit, best, x))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "GBM")
  result <- rbind(result, cc)
  ml.res[["GBM"]] <- list(fit = fit, best = best)
  riskscore[["GBM"]] <- rs

  # survivalsvm
  message("--- 8. survivalsvm ---")
  fit <- train_survivalsvm(est_dd, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "survival-SVM")
  result <- rbind(result, cc)
  ml.res[["survival-SVM"]] <- fit
  riskscore[["survival-SVM"]] <- rs

  # Ridge
  message("--- 9. Ridge ---")
  fit <- train_ridge(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "Ridge")
  result <- rbind(result, cc)
  ml.res[["Ridge"]] <- fit
  riskscore[["Ridge"]] <- rs

  # Lasso
  message("--- 10. Lasso ---")
  fit <- train_lasso(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, "Lasso")
  result <- rbind(result, cc)
  ml.res[["Lasso"]] <- fit
  riskscore[["Lasso"]] <- rs

  # Lasso combinations
  lasso_combos <- c("CoxBoost", "GBM", "plsRcox", "RSF", "StepCox", "SuperPC", "survivalsvm")
  for (algo in lasso_combos) {
    combo_result <- run_lasso_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      pre_var = pre_var,
      seed = seed,
      second_algo = algo
    )
    if (!is.null(combo_result)) {
      result <- rbind(result, combo_result$result)
      ml.res <- c(ml.res, combo_result$ml.res)
      riskscore <- c(riskscore, combo_result$riskscore)
    }
  }

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}

#' Run single algorithm
#'
#' @inheritParams ML.Dev.Prog.Sig
#' @return Result list
#' @keywords internal
run_single_algorithm <- function(est_dd,
                                 train_data,
                                 val_dd_list,
                                 list_train_vali_Data,
                                 pre_var,
                                 single_ml,
                                 rf_nodesize,
                                 seed,
                                 alpha_for_enet,
                                 direction_for_stepcox,
                                 cores_for_parallel) {
  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  if (single_ml == "RSF") {
    message("--- RSF ---")
    set.seed(seed)
    fit <- train_rsf(est_dd, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_rsf(fit, x))
  } else if (single_ml == "Enet") {
    message(paste0("--- Enet[α=", alpha_for_enet, "] ---"))
    set.seed(seed)
    fit <- train_enet(est_dd, pre_var, alpha_for_enet, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
  } else if (single_ml == "StepCox") {
    message(paste0("--- StepCox[", direction_for_stepcox, "] ---"))
    fit <- train_stepcox(est_dd, direction_for_stepcox)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(fit, x))
  } else if (single_ml == "CoxBoost") {
    message("--- CoxBoost ---")
    set.seed(seed)
    fit <- train_coxboost(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(fit, x))
  } else if (single_ml == "plsRcox") {
    message("--- plsRcox ---")
    fit <- train_plsrcox(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
  } else if (single_ml == "SuperPC") {
    message("--- SuperPC ---")
    superpc_result <- train_superpc(est_dd, seed)
    fit <- superpc_result$fit
    cv_fit <- superpc_result$cv_fit
    rs <- lapply(val_dd_list, function(x) {
      cbind(x[, 1:2], RS = predict_superpc(fit, cv_fit, est_dd, x))
    })
  } else if (single_ml == "GBM") {
    message("--- GBM ---")
    gbm_result <- train_gbm(est_dd, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(fit, best, x))
  } else if (single_ml == "survivalsvm") {
    message("--- survivalsvm ---")
    fit <- train_survivalsvm(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
  } else if (single_ml == "Ridge") {
    message("--- Ridge ---")
    fit <- train_ridge(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
  } else if (single_ml == "Lasso") {
    message("--- Lasso ---")
    fit <- train_lasso(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)

  # 构建正确的模型名称（包含参数信息）
  if (single_ml == "StepCox") {
    model_name <- paste0("StepCox[", direction_for_stepcox, "]")
  } else if (single_ml == "Enet") {
    model_name <- paste0("Enet[α=", alpha_for_enet, "]")
  } else {
    model_name <- single_ml
  }

  if (single_ml == "SuperPC") {
    ml.res[[model_name]] <- list(fit, cv_fit)
  } else if (single_ml == "GBM") {
    ml.res[[model_name]] <- list(fit = fit, best = best)
  } else {
    ml.res[[model_name]] <- fit
  }

  riskscore[[model_name]] <- rs
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}

#' Run double algorithm combination
#'
#' @inheritParams ML.Dev.Prog.Sig
#' @return Result list
#' @keywords internal
run_double_algorithm <- function(est_dd,
                                 train_data,
                                 val_dd_list,
                                 list_train_vali_Data,
                                 pre_var,
                                 double_ml1,
                                 double_ml2,
                                 rf_nodesize,
                                 seed,
                                 alpha_for_enet,
                                 direction_for_stepcox,
                                 cores_for_parallel) {
  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Validate algorithm combinations
  if (!double_ml1 %in% c("RSF", "StepCox", "CoxBoost", "Lasso")) {
    stop("double_ml1 must be one of: RSF, StepCox, CoxBoost, Lasso")
  }

  if (!double_ml2 %in% c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", "survivalsvm", "Ridge", "Lasso")) {
    stop("double_ml2 must be a valid algorithm name")
  }

  if (double_ml1 == "RSF") {
    combo_result <- run_rsf_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      rf_nodesize = rf_nodesize,
      seed = seed,
      second_algo = double_ml2,
      alpha_for_enet = alpha_for_enet,
      direction_for_stepcox = direction_for_stepcox,
      cores_for_parallel = cores_for_parallel
    )
  } else if (double_ml1 == "StepCox") {
    combo_result <- run_stepcox_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      direction = direction_for_stepcox,
      seed = seed,
      second_algo = double_ml2,
      cores_for_parallel = cores_for_parallel
    )
  } else if (double_ml1 == "Lasso") {
    combo_result <- run_lasso_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
      pre_var = pre_var,
      seed = seed,
      second_algo = double_ml2,
      direction_for_stepcox = direction_for_stepcox
    )
  }

  if (!is.null(combo_result)) {
    result <- combo_result$result
    ml.res <- combo_result$ml.res
    riskscore <- combo_result$riskscore
  }

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}

# ---- Backward compatibility alias ----

#' @rdname ML.Dev.Prog.Sig
#' @export
ml_dev_prog_sig <- ML.Dev.Prog.Sig
