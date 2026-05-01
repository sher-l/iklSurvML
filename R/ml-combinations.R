# Algorithm Combination Functions for Prognosis
#
# This file contains functions to run combinations of ML algorithms

#' Run RSF + second algorithm combination
#'
#' @param est_dd Training data
#' @param train_data Full training data
#' @param val_dd_list Validation data list
#' @param list_train_vali_Data Original data list for ID recovery
#' @param rf_nodesize Node size for RSF
#' @param seed Random seed
#' @param second_algo Second algorithm name
#' @param alpha_for_enet Alpha for Enet (if applicable)
#' @param direction_for_stepcox Direction for StepCox (if applicable)
#' @param cores_for_parallel Cores for parallel processing
#' @return List with result, ml.res, riskscore
#' @keywords internal
run_rsf_combination <- function(est_dd,
                                train_data,
                                val_dd_list,
                                list_train_vali_Data,
                                rf_nodesize = 5,
                                seed = 5201314,
                                second_algo,
                                alpha_for_enet = NULL,
                                direction_for_stepcox = NULL,
                                cores_for_parallel = 6) {
  valid_algos <- setdiff(survival_ml_names(), "RSF")
  if (!second_algo %in% valid_algos) {
    stop(paste0("RSF combination: second_algo must be one of: ",
                paste(valid_algos, collapse = ", "), ". Got: '", second_algo, "'"))
  }


  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Train RSF and get selected variables
  fit_rsf <- train_rsf(est_dd, rf_nodesize, seed)
  rid <- get_rsf_selected_vars(fit_rsf, seed)

  if (length(rid) <= 1) {
    warning("The number of selected candidate genes by RSF is less than 2")
    return(NULL)
  }

  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", rid)]
  })

  if (second_algo == "CoxBoost") {
    fit <- train_coxboost(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_coxboost(fit, x))
    model_name <- "RSF + CoxBoost"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Enet") {
    alpha_values <- if (is.null(alpha_for_enet)) all_mode_alpha_values() else alpha_for_enet
    for (alpha in alpha_values) {
      fit <- train_enet(est_dd2, rid, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit, x, rid))
      model_name <- paste0("RSF + Enet[\u03b1=", alpha, "]")
      ml.res[[model_name]] <- fit
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      riskscore[[model_name]] <- rs
      cc <- calculate_cindex_result(rs, model_name)
      result <- rbind(result, cc)
    }
    return(list(result = result, ml.res = ml.res, riskscore = riskscore))
  } else if (second_algo == "GBM") {
    gbm_result <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(fit, best, x))
    model_name <- "RSF + GBM"
    ml.res[[model_name]] <- gbm_result
  } else if (second_algo == "Lasso") {
    fit <- train_lasso(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_lasso(fit, x, rid))
    model_name <- "RSF + Lasso"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "plsRcox") {
    fit <- train_plsrcox(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit, x))
    model_name <- "RSF + plsRcox"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Ridge") {
    fit <- train_ridge(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit, x, rid))
    model_name <- "RSF + Ridge"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "StepCox") {
    for (direction in c("both", "backward", "forward")) {
      fit <- train_stepcox(est_dd2, direction)
      rs <- calculate_risk_scores(val_dd_list2, function(x) predict_stepcox(fit, x))
      model_name <- paste0("RSF + StepCox[", direction, "]")
      ml.res[[model_name]] <- fit
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      riskscore[[model_name]] <- rs
      cc <- calculate_cindex_result(rs, model_name)
      result <- rbind(result, cc)
    }
    return(list(result = result, ml.res = ml.res, riskscore = riskscore))
  } else if (second_algo == "SuperPC") {
    superpc_result <- train_superpc(est_dd2, seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(superpc_result, est_dd2, x))
    })
    model_name <- "RSF + SuperPC"
    ml.res[[model_name]] <- superpc_result
  } else if (second_algo == "survivalsvm") {
    fit <- train_survivalsvm(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit, x))
    model_name <- "RSF + survival-SVM"
    ml.res[[model_name]] <- fit
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)
  riskscore[[model_name]] <- rs
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)

  return(list(result = result, ml.res = ml.res, riskscore = riskscore))
}

#' Run StepCox + second algorithm combination
#'
#' @param est_dd Training data
#' @param train_data Full training data
#' @param val_dd_list Validation data list
#' @param list_train_vali_Data Original data list for ID recovery
#' @param direction Direction for StepCox
#' @param seed Random seed
#' @param second_algo Second algorithm name
#' @param alpha_for_enet Alpha for Enet (if applicable)
#' @param cores_for_parallel Cores for parallel processing
#' @return List with result, ml.res, riskscore
#' @keywords internal
run_stepcox_combination <- function(est_dd,
                                    train_data,
                                    val_dd_list,
                                    list_train_vali_Data,
	                                    direction = "both",
	                                    seed = 5201314,
	                                    second_algo,
	                                    alpha_for_enet = NULL,
	                                    rf_nodesize = 5,
	                                    cores_for_parallel = 6) {
  valid_algos <- setdiff(survival_ml_names(), "StepCox")
  if (!second_algo %in% valid_algos) {
    stop(paste0("StepCox combination: second_algo must be one of: ",
                paste(valid_algos, collapse = ", "), ". Got: '", second_algo, "'"))
  }


  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Train StepCox and get selected variables
  fit_stepcox <- train_stepcox(est_dd, direction)
  rid <- get_stepcox_selected_vars(fit_stepcox)

  if (length(rid) <= 1) {
    warning("The number of selected candidate genes by StepCox is less than 2")
    return(NULL)
  }

  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", rid)]
  })

  if (second_algo == "CoxBoost") {
    fit <- train_coxboost(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_coxboost(fit, x))
    model_name <- paste0("StepCox[", direction, "] + CoxBoost")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Enet") {
    alpha_values <- if (is.null(alpha_for_enet)) all_mode_alpha_values() else alpha_for_enet
    for (alpha in alpha_values) {
      fit <- train_enet(est_dd2, rid, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit, x, rid))
      model_name <- paste0("StepCox[", direction, "] + Enet[\u03b1=", alpha, "]")
      ml.res[[model_name]] <- fit
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      riskscore[[model_name]] <- rs
      cc <- calculate_cindex_result(rs, model_name)
      result <- rbind(result, cc)
    }
    return(list(result = result, ml.res = ml.res, riskscore = riskscore))
  } else if (second_algo == "GBM") {
    gbm_result <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(fit, best, x))
    model_name <- paste0("StepCox[", direction, "] + GBM")
    ml.res[[model_name]] <- gbm_result
  } else if (second_algo == "Lasso") {
    fit <- train_lasso(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_lasso(fit, x, rid))
    model_name <- paste0("StepCox[", direction, "] + Lasso")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "plsRcox") {
    fit <- train_plsrcox(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit, x))
    model_name <- paste0("StepCox[", direction, "] + plsRcox")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Ridge") {
    fit <- train_ridge(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit, x, rid))
    model_name <- paste0("StepCox[", direction, "] + Ridge")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "RSF") {
    fit <- train_rsf(est_dd2, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_rsf(fit, x))
    model_name <- paste0("StepCox[", direction, "] + RSF")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "SuperPC") {
    superpc_result <- train_superpc(est_dd2, seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(superpc_result, est_dd2, x))
    })
    model_name <- paste0("StepCox[", direction, "] + SuperPC")
    ml.res[[model_name]] <- superpc_result
  } else if (second_algo == "survivalsvm") {
    fit <- train_survivalsvm(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit, x))
    model_name <- paste0("StepCox[", direction, "] + survival-SVM")
    ml.res[[model_name]] <- fit
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)
  riskscore[[model_name]] <- rs
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)

  return(list(result = result, ml.res = ml.res, riskscore = riskscore))
}

#' Run CoxBoost + second algorithm combination
#'
#' @inheritParams run_rsf_combination
#' @return List with result, ml.res, riskscore
#' @keywords internal
run_coxboost_combination <- function(est_dd,
                                     train_data,
                                     val_dd_list,
                                     list_train_vali_Data,
	                                     seed = 5201314,
	                                     second_algo,
	                                     alpha_for_enet = NULL,
	                                     direction_for_stepcox = "both",
	                                     rf_nodesize = 5,
	                                     cores_for_parallel = 6) {
  valid_algos <- setdiff(survival_ml_names(), "CoxBoost")
  if (!second_algo %in% valid_algos) {
    stop(paste0("CoxBoost combination: second_algo must be one of: ",
                paste(valid_algos, collapse = ", "), ". Got: '", second_algo, "'"))
  }


  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Train CoxBoost and use non-zero coefficients as selected features.
  fit_coxboost <- train_coxboost(est_dd, seed)
  rid <- get_coxboost_selected_vars(fit_coxboost)

  if (length(rid) <= 1) {
    warning("The number of selected candidate genes by CoxBoost is less than 2")
    return(NULL)
  }

  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", rid)]
  })

  if (second_algo == "Enet") {
    alpha_values <- if (is.null(alpha_for_enet)) all_mode_alpha_values() else alpha_for_enet
    for (alpha in alpha_values) {
      fit <- train_enet(est_dd2, rid, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit, x, rid))
      model_name <- paste0("CoxBoost + Enet[\u03b1=", alpha, "]")
      ml.res[[model_name]] <- fit
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      riskscore[[model_name]] <- rs
      cc <- calculate_cindex_result(rs, model_name)
      result <- rbind(result, cc)
    }
    return(list(result = result, ml.res = ml.res, riskscore = riskscore))
  } else if (second_algo == "GBM") {
    gbm_result <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(fit, best, x))
    model_name <- "CoxBoost + GBM"
    ml.res[[model_name]] <- gbm_result
  } else if (second_algo == "Lasso") {
    fit <- train_lasso(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_lasso(fit, x, rid))
    model_name <- "CoxBoost + Lasso"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "plsRcox") {
    fit <- train_plsrcox(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit, x))
    model_name <- "CoxBoost + plsRcox"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Ridge") {
    fit <- train_ridge(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit, x, rid))
    model_name <- "CoxBoost + Ridge"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "RSF") {
    fit <- train_rsf(est_dd2, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_rsf(fit, x))
    model_name <- "CoxBoost + RSF"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "StepCox") {
    fit <- train_stepcox(est_dd2, direction_for_stepcox)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_stepcox(fit, x))
    model_name <- paste0("CoxBoost + StepCox[", direction_for_stepcox, "]")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "SuperPC") {
    superpc_result <- train_superpc(est_dd2, seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(superpc_result, est_dd2, x))
    })
    model_name <- "CoxBoost + SuperPC"
    ml.res[[model_name]] <- superpc_result
  } else if (second_algo == "survivalsvm") {
    fit <- train_survivalsvm(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit, x))
    model_name <- "CoxBoost + survival-SVM"
    ml.res[[model_name]] <- fit
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)
  riskscore[[model_name]] <- rs
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)

  return(list(result = result, ml.res = ml.res, riskscore = riskscore))
}

#' Run Lasso + second algorithm combination
#'
#' @inheritParams run_rsf_combination
#' @return List with result, ml.res, riskscore
#' @keywords internal
run_lasso_combination <- function(est_dd,
                                  train_data,
                                  val_dd_list,
                                  list_train_vali_Data,
	                                  pre_var,
	                                  seed = 5201314,
	                                  second_algo,
	                                  alpha_for_enet = NULL,
	                                  direction_for_stepcox = "both",
	                                  rf_nodesize = 5,
	                                  cores_for_parallel = 6) {
  valid_algos <- setdiff(survival_ml_names(), "Lasso")
  if (!second_algo %in% valid_algos) {
    stop(paste0("Lasso combination: second_algo must be one of: ",
                paste(valid_algos, collapse = ", "), ". Got: '", second_algo, "'"))
  }


  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Train Lasso and get selected variables
  x1 <- as.matrix(est_dd[, pre_var])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  fit_lasso <- glmnet::cv.glmnet(
    x1, x2,
    nfolds = resolve_survival_cv_folds(
      est_dd$OS,
      requested_folds = 10L,
      context = "Lasso selector cross-validation"
    ),
    family = "cox",
    alpha = 1
  )

  rid <- get_lasso_selected_vars(fit_lasso)

  if (length(rid) <= 1) {
    warning("The number of selected candidate genes by Lasso is less than 2")
    return(NULL)
  }

  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", rid)]
  })

  if (second_algo == "CoxBoost") {
    fit <- train_coxboost(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_coxboost(fit, x))
    model_name <- "Lasso + CoxBoost"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "GBM") {
    gbm_result <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(fit, best, x))
    model_name <- "Lasso + GBM"
    ml.res[[model_name]] <- gbm_result
  } else if (second_algo == "plsRcox") {
    fit <- train_plsrcox(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit, x))
    model_name <- "Lasso + plsRcox"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "RSF") {
    fit <- train_rsf(est_dd2, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_rsf(fit, x))
    model_name <- "Lasso + RSF"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "StepCox") {
    fit <- train_stepcox(est_dd2, direction_for_stepcox)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_stepcox(fit, x))
    model_name <- paste0("Lasso + StepCox[", direction_for_stepcox, "]")
    ml.res[[model_name]] <- fit
  } else if (second_algo == "SuperPC") {
    superpc_result <- train_superpc(est_dd2, seed)
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(superpc_result, est_dd2, x))
    })
    model_name <- "Lasso + SuperPC"
    ml.res[[model_name]] <- superpc_result
  } else if (second_algo == "survivalsvm") {
    fit <- train_survivalsvm(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit, x))
    model_name <- "Lasso + survival-SVM"
    ml.res[[model_name]] <- fit
  } else if (second_algo == "Enet") {
    alpha_values <- if (is.null(alpha_for_enet)) all_mode_alpha_values() else alpha_for_enet
    for (alpha in alpha_values) {
      fit <- train_enet(est_dd2, rid, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit, x, rid))
      model_name <- paste0("Lasso + Enet[\u03b1=", alpha, "]")
      ml.res[[model_name]] <- fit
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      riskscore[[model_name]] <- rs
      cc <- calculate_cindex_result(rs, model_name)
      result <- rbind(result, cc)
    }
    return(list(result = result, ml.res = ml.res, riskscore = riskscore))
  } else if (second_algo == "Ridge") {
    fit <- train_ridge(est_dd2, rid, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit, x, rid))
    model_name <- "Lasso + Ridge"
    ml.res[[model_name]] <- fit
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)
  riskscore[[model_name]] <- rs
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)

  return(list(result = result, ml.res = ml.res, riskscore = riskscore))
}

#' Calculate C-index result data frame
#'
#' @param rs Risk score list
#' @param model_name Model name
#' @return Data frame with ID, Cindex, Model columns
#' @keywords internal
calculate_cindex_result <- function(rs, model_name) {
  cc <- data.frame(
    Cindex = sapply(rs, function(x) {
      calculate_cindex(x$RS, x)
    })
  )
  cc$ID <- rownames(cc)
  cc <- cc[, c("ID", "Cindex")]
  cc$Model <- model_name
  return(cc)
}
