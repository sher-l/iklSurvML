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
#'   "SuperPC", "GBM", "survivalsvm", "Ridge", "Lasso".
#' @param alpha_for_Enet Alpha parameter for Enet (0.1 to 0.9).
#' @param direction_for_stepcox Direction for StepCox: "both", "backward", or "forward".
#' @param double_ml1 First algorithm for combination.
#'   Must be one of: "RSF", "StepCox", "CoxBoost", "Lasso".
#' @param double_ml2 Second algorithm for combination.
#'   Must be one of: "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox",
#'   "SuperPC", "GBM", "survivalsvm", "Ridge", "Lasso".
#'   Self-combinations are not supported; choose a different algorithm from
#'   double_ml1 for the second-stage model.
#' @param nodesize Node size parameter for RSF. Default is 5. Try 5-10.
#' @param seed Random seed for reproducibility.
#' @param cores_for_parallel Number of cores for parallel processing. Default is 6.
#' @param feature_alignment How candidate genes are aligned across cohorts.
#'   "strict" (default) requires all candidate genes to be present in every
#'   training/validation cohort. "intersection" preserves the legacy behavior
#'   of silently training only on genes shared by all cohorts, but now warns
#'   about dropped genes.
#' @param allow_partial Logical. For `mode = "all"`, whether to return a
#'   partial model grid when some first-stage selectors or backend models fail.
#'   Defaults to `FALSE` so incomplete all-mode runs fail explicitly instead of
#'   looking complete.
#' @param model_grid All-mode model grid to run. Only `"117"` is supported;
#'   retained for API compatibility and used only when `mode = "all"`.
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
	                            cores_for_parallel = NULL,
	                            feature_alignment = c("strict", "intersection"),
	                            model_grid = "117",
	                            allow_partial = FALSE) {

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

  # ---- Parameter validation ----
  if (is.null(seed)) {
    seed <- 5201314
    message("Using default seed: 5201314")
  }

  if (is.null(rf_nodesize)) {
    rf_nodesize <- 5
  }

  params <- validate_survival_ml_params(
    mode = mode,
    single_ml = single_ml,
    double_ml1 = double_ml1,
    double_ml2 = double_ml2,
    alpha_for_enet = alpha_for_Enet,
    direction_for_stepcox = direction_for_stepcox
  )
  single_ml <- params$single_ml
  double_ml1 <- params$double_ml1
  double_ml2 <- params$double_ml2

	  feature_alignment <- match.arg(feature_alignment)
	  model_grid <- normalize_all_mode_model_grid(model_grid)
	  all_mode_expected <- all_mode_model_grid_size(model_grid)

	  # ---- Data preprocessing ----

  # Normalize '-'/'_' with '.' in feature names
  data_names <- names(list_train_vali_Data)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(list_train_vali_Data))
  }
  list_train_vali_Data <- lapply(data_names, function(nm) {
    normalize_ml_data_columns(list_train_vali_Data[[nm]], paste0("Dataset '", nm, "'"))
  })
  names(list_train_vali_Data) <- data_names

  candidate_genes <- normalize_ml_feature_names(candidate_genes)
  train_data <- normalize_ml_data_columns(train_data, "Training data")

	  # Find/validate features across all datasets without letting validation
	  # cohorts silently shrink the training feature space unless requested.
	  common_feature <- resolve_survival_common_features(
	    train_data = train_data,
	    list_train_vali_Data = list_train_vali_Data,
	    candidate_genes = candidate_genes,
	    feature_alignment = feature_alignment
	  )

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

    # Learn preprocessing from training data only, then apply it to validation cohorts
    preprocessed_train <- preprocess_train_data(train_data, common_feature, return_recipe = TRUE)
    train_data <- preprocessed_train$data
    preprocess_recipe <- preprocessed_train$recipe
    list_train_vali_Data <- preprocess_data_list(
      list_train_vali_Data,
      common_feature,
      recipe = preprocess_recipe
    )

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

      if (length(cd.gene1) < 2) {
        stop(paste0("Insufficient features after univariate Cox filtering: ",
                    length(cd.gene1), ". Need at least 2 features for ML analysis."))
      }

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

    if (length(pre_var) < 2) {
      stop(paste0("Insufficient features after filtering: ", length(pre_var),
                  ". Need at least 2 features for ML analysis."))
    }

    message(paste0("--- Running with ", length(pre_var), " features, ",
                   nrow(est_dd), " training samples, ",
                   length(val_dd_list), " validation cohorts ---"))
    model_warnings <- character()
    capture_model_warnings <- function(expr) {
      withCallingHandlers(
        expr,
        warning = function(w) {
          model_warnings <<- c(model_warnings, conditionMessage(w))
        }
      )
    }

    # ---- Run ML algorithms based on mode ----
    if (mode == "all") {
      result <- capture_model_warnings(run_all_algorithms(
        est_dd = est_dd,
        train_data = train_data,
        val_dd_list = val_dd_list,
        list_train_vali_Data = list_train_vali_Data,
        pre_var = pre_var,
        rf_nodesize = rf_nodesize,
        seed = seed,
        cores_for_parallel = cores_for_parallel,
        model_grid = model_grid
      ))
    } else if (mode == "single") {
      result <- capture_model_warnings(run_single_algorithm(
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
      ))
	    } else if (mode == "double") {
	      result <- capture_model_warnings(run_double_algorithm(
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
	      ))
	    }

	    if (identical(mode, "all")) {
	      assert_complete_all_mode_result(
	        result,
	        expected = all_mode_expected,
	        allow_partial = allow_partial,
	        context = "ML.Dev.Prog.Sig all-mode"
	      )
	    }
	
	    message("--- The analysis has been completed ---")
    if (length(model_warnings) > 0L) {
      result$Model.warnings <- unique(model_warnings)
    }
    result$Preprocess.recipe <- preprocess_recipe
    return(attach_survival_model_info(result))
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
	                               cores_for_parallel,
	                               model_grid = "117") {
  model_grid <- normalize_all_mode_model_grid(model_grid)
  all_mode_expected <- all_mode_model_grid_size(model_grid)
  stepcox_selector_dirs <- all_mode_stepcox_selector_dirs(model_grid)
  coxboost_combos <- all_mode_coxboost_second_stage_algorithms(model_grid)
  lasso_combos <- all_mode_lasso_second_stage_algorithms(model_grid)

  result <- data.frame()
  ml.res <- list()
  riskscore <- list()
  model_errors <- character()
  model_skips <- character()
  record_model_error <- function(model_name, e) {
    msg <- paste0(model_name, ": ", conditionMessage(e))
    model_errors <<- c(model_errors, msg)
    warning(paste0(model_name, " failed: ", conditionMessage(e)), call. = FALSE)
  }
  record_model_skip <- function(model_name, reason) {
    msg <- paste0(model_name, ": ", reason)
    model_skips <<- c(model_skips, msg)
    warning(paste0(model_name, " skipped: ", reason), call. = FALSE)
  }

  # ---- 1. RSF ----
  message("--- 1-1 RSF ---")
  tryCatch({
    set.seed(seed)
    fit <- train_rsf(est_dd, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_rsf(fit, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "RSF")
    result <- rbind(result, cc)
    ml.res[["RSF"]] <- fit
    riskscore[["RSF"]] <- rs
  }, error = function(e) record_model_error("RSF", e))

  # ---- RSF combinations ----
  rsf_combos <- c("CoxBoost", "Enet", "GBM", "Lasso", "plsRcox", "Ridge", "StepCox", "SuperPC", "survivalsvm")
  for (algo in rsf_combos) {
    message(paste0("--- 1-", match(algo, rsf_combos) + 1, " RSF + ", algo, " ---"))
    tryCatch({
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
      } else {
        record_model_skip(paste0("RSF+", algo), "first-stage selector returned fewer than 2 variables")
      }
    }, error = function(e) record_model_error(paste0("RSF+", algo), e))
  }

  # ---- 2. Enet ----
  message("--- 2. Enet ---")
  for (alpha in seq(0.1, 0.9, 0.1)) {
    tryCatch({
      message(paste0("--- 2. Enet[\u03b1=", alpha, "] ---"))
      set.seed(seed)
      fit <- train_enet(est_dd, pre_var, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      cc <- calculate_cindex_result(rs, paste0("Enet[\u03b1=", alpha, "]"))
      result <- rbind(result, cc)
      ml.res[[paste0("Enet[\u03b1=", alpha, "]")]] <- fit
      riskscore[[paste0("Enet[\u03b1=", alpha, "]")]] <- rs
    }, error = function(e) record_model_error(paste0("Enet[alpha=", alpha, "]"), e))
  }

  # ---- 3. StepCox ----
  message("--- 3. StepCox ---")
  for (direction in c("both", "backward", "forward")) {
    tryCatch({
      message(paste0("--- 3. StepCox[", direction, "] ---"))
      fit <- train_stepcox(est_dd, direction)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(fit, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      cc <- calculate_cindex_result(rs, paste0("StepCox[", direction, "]"))
      result <- rbind(result, cc)
      ml.res[[paste0("StepCox[", direction, "]")]] <- fit
      riskscore[[paste0("StepCox[", direction, "]")]] <- rs
    }, error = function(e) record_model_error(paste0("StepCox[", direction, "]"), e))
  }

  # StepCox combinations. The fixed 117-grid includes forward as a first-stage selector.
  stepcox_combos <- c("CoxBoost", "Enet", "GBM", "Lasso", "plsRcox", "Ridge", "RSF", "SuperPC", "survivalsvm")
  for (direction in stepcox_selector_dirs) {
    message(paste0("--- StepCox[", direction, "] combinations ---"))
    for (algo in stepcox_combos) {
      tryCatch({
        combo_result <- run_stepcox_combination(
          est_dd = est_dd,
          train_data = train_data,
          val_dd_list = val_dd_list,
          list_train_vali_Data = list_train_vali_Data,
	          direction = direction,
	          seed = seed,
	          second_algo = algo,
	          rf_nodesize = rf_nodesize,
	          cores_for_parallel = cores_for_parallel
	        )
        if (!is.null(combo_result)) {
          result <- rbind(result, combo_result$result)
          ml.res <- c(ml.res, combo_result$ml.res)
          riskscore <- c(riskscore, combo_result$riskscore)
        } else {
          record_model_skip(paste0("StepCox[", direction, "]+", algo), "first-stage selector returned fewer than 2 variables")
        }
      }, error = function(e) record_model_error(paste0("StepCox[", direction, "]+", algo), e))
    }
  }

  # ---- 4-10. Other single algorithms ----
  # CoxBoost
  message("--- 4. CoxBoost ---")
  tryCatch({
    set.seed(seed)
    fit <- train_coxboost(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(fit, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "CoxBoost")
    result <- rbind(result, cc)
    ml.res[["CoxBoost"]] <- fit
    riskscore[["CoxBoost"]] <- rs
  }, error = function(e) record_model_error("CoxBoost", e))

  # CoxBoost combinations. The fixed 117-grid omits CoxBoost + RSF.
  for (algo in coxboost_combos) {
    if (algo == "StepCox") {
      for (dir in c("both", "backward", "forward")) {
        tryCatch({
          combo_result <- run_coxboost_combination(
            est_dd = est_dd,
            train_data = train_data,
            val_dd_list = val_dd_list,
            list_train_vali_Data = list_train_vali_Data,
	            seed = seed,
	            second_algo = algo,
	            direction_for_stepcox = dir,
	            rf_nodesize = rf_nodesize,
	            cores_for_parallel = cores_for_parallel
	          )
          if (!is.null(combo_result)) {
            result <- rbind(result, combo_result$result)
            ml.res <- c(ml.res, combo_result$ml.res)
            riskscore <- c(riskscore, combo_result$riskscore)
          } else {
            record_model_skip(paste0("CoxBoost+StepCox[", dir, "]"), "first-stage selector returned fewer than 2 variables")
          }
        }, error = function(e) record_model_error(paste0("CoxBoost+StepCox[", dir, "]"), e))
      }
    } else {
      tryCatch({
        combo_result <- run_coxboost_combination(
          est_dd = est_dd,
          train_data = train_data,
          val_dd_list = val_dd_list,
	          list_train_vali_Data = list_train_vali_Data,
	          seed = seed,
	          second_algo = algo,
	          rf_nodesize = rf_nodesize,
	          cores_for_parallel = cores_for_parallel
	        )
        if (!is.null(combo_result)) {
          result <- rbind(result, combo_result$result)
          ml.res <- c(ml.res, combo_result$ml.res)
          riskscore <- c(riskscore, combo_result$riskscore)
        } else {
          record_model_skip(paste0("CoxBoost+", algo), "first-stage selector returned fewer than 2 variables")
        }
      }, error = function(e) record_model_error(paste0("CoxBoost+", algo), e))
    }
  }

  # plsRcox
  message("--- 5. plsRcox ---")
  tryCatch({
    fit <- train_plsrcox(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "plsRcox")
    result <- rbind(result, cc)
    ml.res[["plsRcox"]] <- fit
    riskscore[["plsRcox"]] <- rs
  }, error = function(e) record_model_error("plsRcox", e))

  # SuperPC
  message("--- 6. SuperPC ---")
  tryCatch({
    superpc_result <- train_superpc(est_dd, seed)
    rs <- lapply(val_dd_list, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(superpc_result, est_dd, x))
    })
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "SuperPC")
    result <- rbind(result, cc)
    ml.res[["SuperPC"]] <- superpc_result
    riskscore[["SuperPC"]] <- rs
  }, error = function(e) record_model_error("SuperPC", e))

  # GBM
  message("--- 7. GBM ---")
  tryCatch({
    gbm_result <- train_gbm(est_dd, seed, cores_for_parallel)
    fit <- gbm_result$fit
    best <- gbm_result$best
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(fit, best, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "GBM")
    result <- rbind(result, cc)
    ml.res[["GBM"]] <- list(fit = fit, best = best)
    riskscore[["GBM"]] <- rs
  }, error = function(e) record_model_error("GBM", e))

  # survivalsvm
  message("--- 8. survivalsvm ---")
  tryCatch({
    fit <- train_survivalsvm(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "survival-SVM")
    result <- rbind(result, cc)
    ml.res[["survival-SVM"]] <- fit
    riskscore[["survival-SVM"]] <- rs
  }, error = function(e) record_model_error("survival-SVM", e))

  # Ridge
  message("--- 9. Ridge ---")
  tryCatch({
    fit <- train_ridge(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "Ridge")
    result <- rbind(result, cc)
    ml.res[["Ridge"]] <- fit
    riskscore[["Ridge"]] <- rs
  }, error = function(e) record_model_error("Ridge", e))

  # Lasso
  message("--- 10. Lasso ---")
  tryCatch({
    fit <- train_lasso(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
    rs <- return_id_to_rs(rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, "Lasso")
    result <- rbind(result, cc)
    ml.res[["Lasso"]] <- fit
    riskscore[["Lasso"]] <- rs
  }, error = function(e) record_model_error("Lasso", e))

  # Lasso combinations. The fixed 117-grid omits elastic-net and ridge after Lasso selection.
  for (algo in lasso_combos) {
    if (algo == "StepCox") {
      for (dir in c("both", "backward", "forward")) {
        tryCatch({
          combo_result <- run_lasso_combination(
            est_dd = est_dd,
            train_data = train_data,
            val_dd_list = val_dd_list,
            list_train_vali_Data = list_train_vali_Data,
            pre_var = pre_var,
	            seed = seed,
	            second_algo = algo,
	            direction_for_stepcox = dir,
	            rf_nodesize = rf_nodesize,
	            cores_for_parallel = cores_for_parallel
	          )
          if (!is.null(combo_result)) {
            result <- rbind(result, combo_result$result)
            ml.res <- c(ml.res, combo_result$ml.res)
            riskscore <- c(riskscore, combo_result$riskscore)
          } else {
            record_model_skip(paste0("Lasso+StepCox[", dir, "]"), "first-stage selector returned fewer than 2 variables")
          }
        }, error = function(e) record_model_error(paste0("Lasso+StepCox[", dir, "]"), e))
      }
    } else {
      tryCatch({
        combo_result <- run_lasso_combination(
          est_dd = est_dd,
          train_data = train_data,
          val_dd_list = val_dd_list,
	          list_train_vali_Data = list_train_vali_Data,
	          pre_var = pre_var,
	          seed = seed,
	          second_algo = algo,
	          rf_nodesize = rf_nodesize,
	          cores_for_parallel = cores_for_parallel
	        )
        if (!is.null(combo_result)) {
          result <- rbind(result, combo_result$result)
          ml.res <- c(ml.res, combo_result$ml.res)
          riskscore <- c(riskscore, combo_result$riskscore)
        } else {
          record_model_skip(paste0("Lasso+", algo), "first-stage selector returned fewer than 2 variables")
        }
      }, error = function(e) record_model_error(paste0("Lasso+", algo), e))
    }
  }

  if (length(model_errors) > 0 && length(ml.res) == 0L) {
    stop(paste0(
      "Model training failed for all all-mode tasks: ",
      paste(model_errors, collapse = "; ")
    ), call. = FALSE)
  }
  if (length(model_errors) > 0) {
    warning(paste0(
      "Model training failed for ",
      length(model_errors),
      " all-mode task(s); returning successfully fitted models. Failed tasks: ",
      paste(model_errors, collapse = "; ")
    ), call. = FALSE)
  }

  warn_if_all_mode_incomplete(length(ml.res), expected = all_mode_expected, context = "All-mode")

  out <- list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  )
  if (length(model_skips) > 0L) {
    out$Model.skips <- unique(model_skips)
  }
  if (length(model_errors) > 0L) {
    out$Model.errors <- unique(model_errors)
  }
  return(out)
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
    message(paste0("--- Enet[\u03b1=", alpha_for_enet, "] ---"))
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
    fit <- train_superpc(est_dd, seed)
    rs <- lapply(val_dd_list, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(fit, est_dd, x))
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
    model_name <- paste0("Enet[\u03b1=", alpha_for_enet, "]")
  } else {
    model_name <- display_survival_ml_name(single_ml)
  }

  if (single_ml == "SuperPC") {
    ml.res[[model_name]] <- fit
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
  if (!double_ml1 %in% survival_ml_first_stage_names()) {
    stop("double_ml1 must be one of: RSF, StepCox, CoxBoost, Lasso")
  }

  if (!double_ml2 %in% survival_ml_names()) {
    stop("double_ml2 must be a valid algorithm name")
  }
  assert_no_survival_self_combination(double_ml1, double_ml2)

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
	      rf_nodesize = rf_nodesize,
	      cores_for_parallel = cores_for_parallel
	    )
  } else if (double_ml1 == "CoxBoost") {
    combo_result <- run_coxboost_combination(
      est_dd = est_dd,
      train_data = train_data,
      val_dd_list = val_dd_list,
      list_train_vali_Data = list_train_vali_Data,
	      seed = seed,
	      second_algo = double_ml2,
	      direction_for_stepcox = direction_for_stepcox,
	      rf_nodesize = rf_nodesize,
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
	      direction_for_stepcox = direction_for_stepcox,
	      rf_nodesize = rf_nodesize,
	      cores_for_parallel = cores_for_parallel
	    )
  }

  if (!is.null(combo_result)) {
    result <- combo_result$result
    ml.res <- combo_result$ml.res
    riskscore <- combo_result$riskscore
  } else {
    stop(
      "The first-stage algorithm selected fewer than 2 variables; double-mode model was not fit.",
      call. = FALSE
    )
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
