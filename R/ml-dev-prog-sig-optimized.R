# Optimized Main Function - Fixed All-Mode Combinations
#
# Key optimization: Train feature selectors once and reuse for all combinations
# Fixed Mime-style all-mode grid: 117 model combinations

#' Optimized version of ML.Dev.Prog.Sig with caching
#'
#' @description
#' Same functionality as ML.Dev.Prog.Sig but with significant performance improvements:
#' - Feature selection algorithms (RSF, StepCox, CoxBoost, Lasso) trained once and reused
#' - Fixed 117 all-mode combinations
#' - Optional forked parallel execution for all-mode combinations (Linux/macOS fork)
#'
#' @inheritParams ML.Dev.Prog.Sig
#' @param use_parallel Enable forked parallel processing for mode="all".
#'   Defaults to FALSE because sequential all-mode is the safest cross-platform
#'   path; TRUE keeps GBM tasks in the parent process for fork safety.
#' @param feature_alignment How candidate genes are aligned across cohorts.
#'   "strict" (default) requires all candidate genes to be present in every
#'   training/validation cohort. "intersection" preserves the legacy behavior
#'   of training only on genes shared by all cohorts, with a warning for drops.
#' @param allow_partial Logical. For `mode = "all"`, whether to return a
#'   partial model grid when some first-stage selectors or backend models fail.
#'   Defaults to `FALSE` so incomplete all-mode runs fail explicitly instead of
#'   looking complete.
#' @param fallback_on_fast_failure Logical. For `mode = "all"`, whether to
#'   retry with the resilient legacy all-mode runner after the optimized runner
#'   fails before returning results. Defaults to `FALSE` so optimized failures
#'   are explicit.
#' @return Same result structure as ML.Dev.Prog.Sig
#' @export
ML.Dev.Prog.Sig.Fast <- function(train_data,
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
	                                  cores_for_parallel = 12,
	                                  use_parallel = FALSE,
	                                  feature_alignment = c("strict", "intersection"),
	                                  model_grid = "117",
	                                  allow_partial = FALSE,
	                                  fallback_on_fast_failure = FALSE) {

  # ---- Set default parameters ----
  if (is.null(alpha_for_Enet)) alpha_for_Enet <- 0.1
  if (is.null(direction_for_stepcox)) direction_for_stepcox <- "both"
  if (is.null(unicox_p_cutoff)) unicox_p_cutoff <- 0.05
  if (is.null(unicox.filter.for.candi)) unicox.filter.for.candi <- TRUE
  if (is.null(seed)) {
    seed <- 5201314
    message("Using default seed: 5201314")
  }

  rf_nodesize <- if(is.null(nodesize)) 5 else nodesize
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

	  # ---- Data preprocessing (same as original) ----
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

	  common_feature <- resolve_survival_common_features(
	    train_data = train_data,
	    list_train_vali_Data = list_train_vali_Data,
	    candidate_genes = candidate_genes,
	    feature_alignment = feature_alignment
	  )

  # Validate parameters
  if (length(rf_nodesize) != 1 || is.na(rf_nodesize) ||
      length(seed) != 1 || is.na(seed)) {
    stop("nodesize and seed must be non-missing scalar values", call. = FALSE)
  }
  if (!is.logical(fallback_on_fast_failure) ||
      length(fallback_on_fast_failure) != 1L ||
      is.na(fallback_on_fast_failure)) {
    stop("fallback_on_fast_failure must be TRUE or FALSE", call. = FALSE)
  }
  if (!identical(c("ID", "OS.time", "OS"), colnames(train_data)[1:3])) {
    stop("first 3 columns of train_data must be ID, OS.time, OS", call. = FALSE)
  }
  if (length(candidate_genes) == 0 || length(common_feature) <= 3) {
    stop("candidate_genes must exist across train_data and all validation datasets", call. = FALSE)
  }

  # Learn preprocessing from training data only, then apply it to validation cohorts
  preprocessed_train <- preprocess_train_data(train_data, common_feature, return_recipe = TRUE)
  train_data <- preprocessed_train$data
  preprocess_recipe <- preprocessed_train$recipe
  list_train_vali_Data <- preprocess_data_list(
    list_train_vali_Data,
    common_feature,
    recipe = preprocess_recipe
  )

  # Univariate Cox filtering
  if (unicox.filter.for.candi) {
    cd.gene <- common_feature[-c(1:3)]
    cd.gene1 <- sig_unicox(gene_list = cd.gene, inputSet = train_data, unicox_pcutoff = unicox_p_cutoff)
    message(paste0("--- Unicox filtered genes: ", length(cd.gene1), " ---"))
    if (length(cd.gene1) < 2) {
      stop(paste0("Insufficient features after univariate Cox filtering: ",
                  length(cd.gene1), ". Need at least 2 features for ML analysis."), call. = FALSE)
    }
    common_feature <- c(common_feature[1:3], cd.gene1)
    train_data <- as.data.frame(train_data)[, common_feature]
  }

  est_dd <- as.data.frame(train_data)[, common_feature[-1]]
  val_dd_list <- lapply(list_train_vali_Data, function(x) x[, common_feature[-1]])
  pre_var <- common_feature[-c(1:3)]
  model_warnings <- character()
  capture_model_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        model_warnings <<- c(model_warnings, conditionMessage(w))
      }
    )
  }

  # ---- Run algorithms ----
  if (mode == "all") {
	    result <- capture_model_warnings({
	      fast_runner <- if (use_parallel) {
	        function() run_all_algorithms_128_parallel(
	          est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
	          rf_nodesize, seed, cores = cores_for_parallel, model_grid = model_grid
	        )
	      } else {
	        function() run_all_algorithms_128(
	          est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
	          rf_nodesize, seed, cores_for_parallel, model_grid = model_grid
	        )
	      }

      tryCatch(
        fast_runner(),
        error = function(e) {
          original_error <- conditionMessage(e)
          if (!isTRUE(fallback_on_fast_failure)) {
            stop(paste0(
              "Fast all-mode failed before returning model results: ",
              original_error,
              ". Set fallback_on_fast_failure = TRUE to retry with the resilient all-mode runner."
            ), call. = FALSE)
          }
          warning(paste0(
            "Fast all-mode failed before returning model results; falling back ",
            "to the resilient all-mode runner. Original error: ",
            original_error
          ), call. = FALSE)
          fallback <- run_all_algorithms(
            est_dd = est_dd,
            train_data = train_data,
            val_dd_list = val_dd_list,
            list_train_vali_Data = list_train_vali_Data,
	            pre_var = pre_var,
	            rf_nodesize = rf_nodesize,
	            seed = seed,
	            cores_for_parallel = cores_for_parallel,
	            model_grid = model_grid
	          )
          fallback$Fast.fallback <- list(
            original_error = original_error,
            fallback_runner = "run_all_algorithms",
            optimized_runner = if (use_parallel) {
              "run_all_algorithms_128_parallel"
            } else {
              "run_all_algorithms_128"
            }
          )
          fallback
	        }
	      )
    })
  } else if (mode == "single") {
    result <- capture_model_warnings(run_single_algorithm(
      est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
      single_ml, rf_nodesize, seed, alpha_for_Enet, direction_for_stepcox, cores_for_parallel
    ))
	  } else if (mode == "double") {
	    result <- capture_model_warnings(run_double_algorithm_optimized(
	      est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
	      double_ml1, double_ml2, rf_nodesize, seed, alpha_for_Enet,
	      direction_for_stepcox, cores_for_parallel
	    ))
	  }

	  if (identical(mode, "all")) {
		    assert_complete_all_mode_result(
		      result,
		      expected = all_mode_expected,
		      allow_partial = allow_partial,
		      context = "ML.Dev.Prog.Sig.Fast all-mode"
		    )
	  }
	
	  message("--- Analysis completed ---")
  if (length(model_warnings) > 0L) {
    result$Model.warnings <- unique(model_warnings)
  }
  result$Preprocess.recipe <- preprocess_recipe
  return(attach_survival_model_info(result))
}

#' Helper function to add model result
#' @keywords internal
add_model_result <- function(result, ml.res, riskscore, rs, fit, model_name, list_train_vali_Data) {
  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)
  ml.res[[model_name]] <- fit
  riskscore[[model_name]] <- rs
  list(result = result, ml.res = ml.res, riskscore = riskscore)
}

#' Fixed 117 all-mode combinations - Optimized
#' @keywords internal
run_all_algorithms_128 <- function(est_dd, train_data, val_dd_list,
                                    list_train_vali_Data, pre_var,
                                    rf_nodesize, seed, cores_for_parallel,
                                    model_grid = "117") {
  model_grid <- normalize_all_mode_model_grid(model_grid)
  all_mode_expected <- all_mode_model_grid_size(model_grid)
  stepcox_selector_dirs <- all_mode_stepcox_selector_dirs(model_grid)

  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # ============================================
  # PHASE 1: Pre-compute feature selectors
  # ============================================
  message("--- Phase 1: Training feature selectors (RSF, StepCox, CoxBoost, Lasso) ---")

  # 1. RSF - train once, get selected vars
  rsf_fit <- train_rsf(est_dd, rf_nodesize, seed)
  set.seed(seed)
  rsf_vars <- randomForestSRC::var.select(object = rsf_fit, conservative = "high")$topvars
  message(paste0("  RSF selected ", length(rsf_vars), " variables"))

  # 2. StepCox - train all 3 directions once
  stepcox_fits <- list()
  stepcox_vars <- list()
  for (direction in c("both", "backward", "forward")) {
    stepcox_fits[[direction]] <- train_stepcox(est_dd, direction)
    stepcox_vars[[direction]] <- get_stepcox_selected_vars(stepcox_fits[[direction]])
    message(paste0("  StepCox[", direction, "] selected ", length(stepcox_vars[[direction]]), " variables"))
  }

  # 3. CoxBoost - train once
  coxboost_fit <- train_coxboost(est_dd, seed)
  coxboost_vars <- get_coxboost_selected_vars(coxboost_fit)
  message(paste0("  CoxBoost selected ", length(coxboost_vars), " variables"))

  # 4. Lasso - train once
  lasso_fit <- train_lasso(est_dd, pre_var, seed)
  lasso_vars <- get_lasso_selected_vars(lasso_fit)
  message(paste0("  Lasso selected ", length(lasso_vars), " variables"))

  # ============================================
  # PHASE 2: Single models (20)
  # ============================================
  message("--- Phase 2: Single models (20) ---")

  # RSF (1)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict(rsf_fit, newdata = x)$predicted)
  tmp <- add_model_result(result, ml.res, riskscore, rs, rsf_fit, "RSF", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # Enet (9 alphas)
  for (alpha in all_mode_alpha_values()) {
    fit <- train_enet(est_dd, pre_var, alpha, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("Enet[\u03b1=", alpha, "]"), list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # StepCox (3 directions)
  for (direction in c("both", "backward", "forward")) {
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(stepcox_fits[[direction]], x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, stepcox_fits[[direction]],
                            paste0("StepCox[", direction, "]"), list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # CoxBoost (1)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(coxboost_fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, coxboost_fit, "CoxBoost", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # Other singles (5): plsRcox, SuperPC, GBM, survival-SVM, Ridge
  fit <- train_plsrcox(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "plsRcox", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  superpc_res <- train_superpc(est_dd, seed)
  rs <- lapply(val_dd_list, function(x) cbind(x[, 1:2], RS = predict_superpc_model(superpc_res, est_dd, x)))
  tmp <- add_model_result(result, ml.res, riskscore, rs, superpc_res, "SuperPC", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  gbm_res <- train_gbm(est_dd, seed, cores_for_parallel)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm_res$fit, gbm_res$best, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_res, "GBM", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  fit <- train_survivalsvm(est_dd, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "survival-SVM", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  fit <- train_ridge(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Ridge", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # Lasso (1)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(lasso_fit, x, pre_var))
  tmp <- add_model_result(result, ml.res, riskscore, rs, lasso_fit, "Lasso", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # ============================================
  # PHASE 3: RSF combinations (19)
  # ============================================
  message("--- Phase 3: RSF combinations (19) ---")

  if (length(rsf_vars) > 1) {
    est_dd_rsf <- train_data[, c("OS.time", "OS", rsf_vars)]
    val_dd_list_rsf <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", rsf_vars)])

    # RSF + CoxBoost
    fit <- train_coxboost(est_dd_rsf, seed)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_coxboost(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "RSF + CoxBoost", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # RSF + Enet (9)
    for (alpha in all_mode_alpha_values()) {
      fit <- train_enet(est_dd_rsf, rsf_vars, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_enet(fit, x, rsf_vars))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("RSF + Enet[\u03b1=", alpha, "]"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }

    # RSF + GBM, Lasso, plsRcox, Ridge, SuperPC, survival-SVM (6)
    gbm_r <- train_gbm(est_dd_rsf, seed, cores_for_parallel)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_gbm(gbm_r$fit, gbm_r$best, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_r, "RSF + GBM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    fit <- train_lasso(est_dd_rsf, rsf_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_lasso(fit, x, rsf_vars))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "RSF + Lasso", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    fit <- train_plsrcox(est_dd_rsf, rsf_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_plsrcox(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "RSF + plsRcox", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    fit <- train_ridge(est_dd_rsf, rsf_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_ridge(fit, x, rsf_vars))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "RSF + Ridge", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    spc_r <- train_superpc(est_dd_rsf, seed)
    rs <- lapply(val_dd_list_rsf, function(x) cbind(x[, 1:2], RS = predict_superpc_model(spc_r, est_dd_rsf, x)))
    tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "RSF + SuperPC", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    fit <- train_survivalsvm(est_dd_rsf, seed)
    rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_survivalsvm(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "RSF + survival-SVM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # RSF + StepCox (3 directions)
    for (direction in c("both", "backward", "forward")) {
      fit <- train_stepcox(est_dd_rsf, direction)
      rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_stepcox(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("RSF + StepCox[", direction, "]"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }
  }

  # ============================================
  # PHASE 4: StepCox combinations (117-grid includes forward as first-stage selector)
  # ============================================
  message(paste0("--- Phase 4: StepCox combinations (", length(stepcox_selector_dirs) * 17L, ") ---"))

  for (direction in stepcox_selector_dirs) {
    sc_vars <- stepcox_vars[[direction]]
    if (length(sc_vars) > 1) {
      est_dd_sc <- train_data[, c("OS.time", "OS", sc_vars)]
      val_dd_list_sc <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", sc_vars)])

      prefix <- paste0("StepCox[", direction, "]")

      # StepCox + CoxBoost
      fit <- train_coxboost(est_dd_sc, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_coxboost(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + CoxBoost"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + Enet (9)
      for (alpha in all_mode_alpha_values()) {
        fit <- train_enet(est_dd_sc, sc_vars, alpha, seed)
        rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_enet(fit, x, sc_vars))
        tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + Enet[\u03b1=", alpha, "]"), list_train_vali_Data)
        result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
      }

      # StepCox + GBM
      gbm_r <- train_gbm(est_dd_sc, seed, cores_for_parallel)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_gbm(gbm_r$fit, gbm_r$best, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_r, paste0(prefix, " + GBM"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + Lasso
      fit <- train_lasso(est_dd_sc, sc_vars, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_lasso(fit, x, sc_vars))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + Lasso"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + plsRcox
      fit <- train_plsrcox(est_dd_sc, sc_vars, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_plsrcox(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + plsRcox"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + Ridge
      fit <- train_ridge(est_dd_sc, sc_vars, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_ridge(fit, x, sc_vars))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + Ridge"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + RSF
      fit <- train_rsf(est_dd_sc, rf_nodesize, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_rsf(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + RSF"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + SuperPC
      spc_r <- train_superpc(est_dd_sc, seed)
      rs <- lapply(val_dd_list_sc, function(x) cbind(x[, 1:2], RS = predict_superpc_model(spc_r, est_dd_sc, x)))
      tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, paste0(prefix, " + SuperPC"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + survival-SVM
      fit <- train_survivalsvm(est_dd_sc, seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_survivalsvm(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + survival-SVM"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }
  }

  # ============================================
  # PHASE 5: CoxBoost combinations (18; 117-grid omits CoxBoost + RSF)
  # ============================================
  message("--- Phase 5: CoxBoost combinations (18) ---")

  if (length(coxboost_vars) > 1) {
    est_dd_coxboost <- train_data[, c("OS.time", "OS", coxboost_vars)]
    val_dd_list_coxboost <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", coxboost_vars)])

    # CoxBoost + Enet (9)
    for (alpha in all_mode_alpha_values()) {
      fit <- train_enet(est_dd_coxboost, coxboost_vars, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_enet(fit, x, coxboost_vars))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("CoxBoost + Enet[\u03b1=", alpha, "]"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }

    # CoxBoost + GBM
    gbm_r <- train_gbm(est_dd_coxboost, seed, cores_for_parallel)
    rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_gbm(gbm_r$fit, gbm_r$best, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_r, "CoxBoost + GBM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # CoxBoost + Lasso
    fit <- train_lasso(est_dd_coxboost, coxboost_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_lasso(fit, x, coxboost_vars))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + Lasso", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # CoxBoost + plsRcox
    fit <- train_plsrcox(est_dd_coxboost, coxboost_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_plsrcox(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + plsRcox", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # CoxBoost + Ridge
    fit <- train_ridge(est_dd_coxboost, coxboost_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_ridge(fit, x, coxboost_vars))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + Ridge", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # CoxBoost + StepCox (3 directions)
    for (direction in c("both", "backward", "forward")) {
      fit <- train_stepcox(est_dd_coxboost, direction)
      rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_stepcox(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("CoxBoost + StepCox[", direction, "]"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }

    # CoxBoost + SuperPC
    spc_r <- train_superpc(est_dd_coxboost, seed)
    rs <- lapply(val_dd_list_coxboost, function(x) cbind(x[, 1:2], RS = predict_superpc_model(spc_r, est_dd_coxboost, x)))
    tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "CoxBoost + SuperPC", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # CoxBoost + survival-SVM
    fit <- train_survivalsvm(est_dd_coxboost, seed)
    rs <- calculate_risk_scores(val_dd_list_coxboost, function(x) predict_survivalsvm(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + survival-SVM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  } else {
    warning("The number of selected candidate genes by CoxBoost is less than 2")
  }

  # ============================================
  # PHASE 6: Lasso combinations (117-grid omits elastic-net and ridge after Lasso selection)
  # ============================================
  message("--- Phase 6: Lasso combinations (9) ---")

  if (length(lasso_vars) > 1) {
    est_dd_lasso <- train_data[, c("OS.time", "OS", lasso_vars)]
    val_dd_list_lasso <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", lasso_vars)])

    # Lasso + CoxBoost
    fit <- train_coxboost(est_dd_lasso, seed)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_coxboost(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Lasso + CoxBoost", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + GBM
    gbm_r <- train_gbm(est_dd_lasso, seed, cores_for_parallel)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_gbm(gbm_r$fit, gbm_r$best, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_r, "Lasso + GBM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + plsRcox
    fit <- train_plsrcox(est_dd_lasso, lasso_vars, seed)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_plsrcox(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Lasso + plsRcox", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + RSF
    fit <- train_rsf(est_dd_lasso, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_rsf(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Lasso + RSF", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + StepCox (3 directions)
    for (direction in c("both", "backward", "forward")) {
      fit <- train_stepcox(est_dd_lasso, direction)
      rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_stepcox(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("Lasso + StepCox[", direction, "]"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }

    # Lasso + SuperPC
    spc_r <- train_superpc(est_dd_lasso, seed)
    rs <- lapply(val_dd_list_lasso, function(x) cbind(x[, 1:2], RS = predict_superpc_model(spc_r, est_dd_lasso, x)))
    tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "Lasso + SuperPC", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + survival-SVM
    fit <- train_survivalsvm(est_dd_lasso, seed)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_survivalsvm(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Lasso + survival-SVM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # ============================================
  # Summary
  # ============================================
  message(paste0("--- Total models: ", length(ml.res), " ---"))
  warn_if_all_mode_incomplete(length(ml.res), expected = all_mode_expected, context = "All-mode")

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}

#' Optimized double algorithm runner
#' @keywords internal
run_double_algorithm_optimized <- function(est_dd, train_data, val_dd_list,
                                            list_train_vali_Data, pre_var,
                                            double_ml1, double_ml2, rf_nodesize, seed,
                                            alpha_for_enet, direction_for_stepcox, cores_for_parallel) {
  assert_no_survival_self_combination(double_ml1, double_ml2)

  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Get first algorithm's selected variables
  if (double_ml1 == "RSF") {
    fit1 <- train_rsf(est_dd, rf_nodesize, seed)
    set.seed(seed)
    selected_vars <- randomForestSRC::var.select(object = fit1, conservative = "high")$topvars
  } else if (double_ml1 == "StepCox") {
    fit1 <- train_stepcox(est_dd, direction_for_stepcox)
    selected_vars <- get_stepcox_selected_vars(fit1)
    double_ml1_display <- paste0("StepCox[", direction_for_stepcox, "]")
  } else if (double_ml1 == "CoxBoost") {
    fit1 <- train_coxboost(est_dd, seed)
    selected_vars <- get_coxboost_selected_vars(fit1)
  } else if (double_ml1 == "Lasso") {
    fit1 <- train_lasso(est_dd, pre_var, seed)
    selected_vars <- get_lasso_selected_vars(fit1)
  }

  # 构建第一个算法的显示名称（StepCox 需要包含方向信息）
  if (!exists("double_ml1_display")) {
    double_ml1_display <- double_ml1
  }

  if (length(selected_vars) <= 1) {
    stop(
      "The first-stage algorithm selected fewer than 2 variables; double-mode model was not fit.",
      call. = FALSE
    )
  }

  est_dd2 <- train_data[, c("OS.time", "OS", selected_vars)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", selected_vars)]
  })

  # Run second algorithm
  if (double_ml2 == "CoxBoost") {
    fit2 <- train_coxboost(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_coxboost(fit2, x))
    model_name <- paste0(double_ml1_display, " + CoxBoost")
  } else if (double_ml2 == "Enet") {
    fit2 <- train_enet(est_dd2, selected_vars, alpha_for_enet, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit2, x, selected_vars))
    model_name <- paste0(double_ml1_display, " + Enet[\u03b1=", alpha_for_enet, "]")
  } else if (double_ml2 == "GBM") {
    gbm_res <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit2 <- gbm_res
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(gbm_res$fit, gbm_res$best, x))
    model_name <- paste0(double_ml1_display, " + GBM")
  } else if (double_ml2 == "Lasso") {
    fit2 <- train_lasso(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_lasso(fit2, x, selected_vars))
    model_name <- paste0(double_ml1_display, " + Lasso")
  } else if (double_ml2 == "plsRcox") {
    fit2 <- train_plsrcox(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit2, x))
    model_name <- paste0(double_ml1_display, " + plsRcox")
  } else if (double_ml2 == "Ridge") {
    fit2 <- train_ridge(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit2, x, selected_vars))
    model_name <- paste0(double_ml1_display, " + Ridge")
  } else if (double_ml2 == "RSF") {
    fit2 <- train_rsf(est_dd2, rf_nodesize, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_rsf(fit2, x))
    model_name <- paste0(double_ml1_display, " + RSF")
  } else if (double_ml2 == "StepCox") {
    fit2 <- train_stepcox(est_dd2, direction_for_stepcox)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_stepcox(fit2, x))
    model_name <- paste0(double_ml1_display, " + StepCox[", direction_for_stepcox, "]")
  } else if (double_ml2 == "SuperPC") {
    spc_res <- train_superpc(est_dd2, seed)
    fit2 <- spc_res
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc_model(spc_res, est_dd2, x))
    })
    model_name <- paste0(double_ml1_display, " + SuperPC")
  } else if (double_ml2 == "survivalsvm") {
    fit2 <- train_survivalsvm(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit2, x))
    model_name <- paste0(double_ml1_display, " + survival-SVM")
  }

  rs <- return_id_to_rs(rs, list_train_vali_Data)
  cc <- calculate_cindex_result(rs, model_name)
  result <- rbind(result, cc)
  ml.res[[model_name]] <- fit2
  riskscore[[model_name]] <- rs

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}
