# Optimized Main Function - Complete 117 Combinations
#
# Key optimization: Train feature selectors once and reuse for all combinations
# Matches original Mime package exactly: 117 model combinations

#' Optimized version of ML.Dev.Prog.Sig with caching - Complete 117 combinations
#'
#' @description
#' Same functionality as ML.Dev.Prog.Sig but with significant performance improvements:
#' - Feature selection algorithms (RSF, StepCox, CoxBoost, Lasso) trained once and reused
#' - All 117 combinations as per the original Mime paper
#'
#' @inheritParams ML.Dev.Prog.Sig
#' @param use_parallel Enable parallel processing (requires future package)
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
                                  cores_for_parallel = NULL,
                                  use_parallel = TRUE) {

  # ---- Set default parameters ----
  if (is.null(alpha_for_Enet)) alpha_for_Enet <- 0.1
  if (is.null(cores_for_parallel)) cores_for_parallel <- parallelly::availableCores(omit = 1)
  if (is.null(direction_for_stepcox)) direction_for_stepcox <- "both"
  if (is.null(unicox_p_cutoff)) unicox_p_cutoff <- 0.05
  if (is.null(unicox.filter.for.candi)) unicox.filter.for.candi <- TRUE

  rf_nodesize <- if(is.null(nodesize)) 5 else nodesize

  # ---- Setup parallel backend ----
  if (use_parallel && requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = cores_for_parallel)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # ---- Data preprocessing (same as original) ----
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    x
  })
  candidate_genes <- gsub("-", ".", candidate_genes)
  colnames(train_data) <- gsub("-", ".", colnames(train_data))

  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  # Validate parameters
  if (!(!is.na(rf_nodesize) && !is.na(seed) &&
        mode %in% c("all", "single", "double") &&
        identical(c("ID", "OS.time", "OS"), colnames(train_data)[1:3]) &&
        length(candidate_genes) > 0 &&
        length(common_feature) > 3)) {
    stop("Invalid parameters or data format")
  }

  # Preprocess data
  list_train_vali_Data <- preprocess_data_list(list_train_vali_Data, common_feature)
  train_data <- preprocess_train_data(train_data, common_feature)

  # Univariate Cox filtering
  if (unicox.filter.for.candi) {
    cd.gene <- common_feature[-c(1:3)]
    cd.gene1 <- sig_unicox(gene_list = cd.gene, inputSet = train_data, unicox_pcutoff = unicox_p_cutoff)
    message(paste0("--- Unicox filtered genes: ", length(cd.gene1), " ---"))
    common_feature <- c(common_feature[1:3], cd.gene1)
    train_data <- as.data.frame(train_data)[, common_feature]
  }

  est_dd <- as.data.frame(train_data)[, common_feature[-1]]
  val_dd_list <- lapply(list_train_vali_Data, function(x) x[, common_feature[-1]])
  pre_var <- common_feature[-c(1:3)]

  # ---- Run algorithms ----
  if (mode == "all") {
    result <- run_all_algorithms_117(
      est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
      rf_nodesize, seed, cores_for_parallel
    )
  } else if (mode == "single") {
    result <- run_single_algorithm(
      est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
      single_ml, rf_nodesize, seed, alpha_for_Enet, direction_for_stepcox, cores_for_parallel
    )
  } else if (mode == "double") {
    result <- run_double_algorithm_optimized(
      est_dd, train_data, val_dd_list, list_train_vali_Data, pre_var,
      double_ml1, double_ml2, rf_nodesize, seed, alpha_for_Enet,
      direction_for_stepcox, cores_for_parallel
    )
  }

  message("--- Analysis completed ---")
  return(result)
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

#' Complete 117 combinations - Optimized
#' @keywords internal
run_all_algorithms_117 <- function(est_dd, train_data, val_dd_list,
                                    list_train_vali_Data, pre_var,
                                    rf_nodesize, seed, cores_for_parallel) {
  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # ============================================
  # PHASE 1: Pre-compute feature selectors
  # ============================================
  message("--- Phase 1: Training feature selectors (RSF, StepCox, CoxBoost, Lasso) ---")

  # 1. RSF - train once, get selected vars
  set.seed(seed)
  Surv <- survival::Surv  # Fix scoping
  rsf_fit <- randomForestSRC::rfsrc(
    Surv(OS.time, OS) ~ .,
    data = est_dd,
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = "logrank",
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE,
    seed = seed
  )
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
  # CoxBoost doesn't have explicit feature selection, use all features
  coxboost_vars <- pre_var

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
  for (alpha in seq(0.1, 0.9, 0.1)) {
    fit <- train_enet(est_dd, pre_var, alpha, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("Enet[α=", alpha, "]"), list_train_vali_Data)
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
  # Note: Original code uses "survival - SVM" for single model (with spaces)
  fit <- train_plsrcox(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "plsRcox", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  superpc_res <- train_superpc(est_dd, seed)
  rs <- lapply(val_dd_list, function(x) cbind(x[, 1:2], RS = predict_superpc(superpc_res$fit, superpc_res$cv_fit, est_dd, x)))
  tmp <- add_model_result(result, ml.res, riskscore, rs, superpc_res, "SuperPC", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  gbm_res <- train_gbm(est_dd, seed, cores_for_parallel)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm_res$fit, gbm_res$best, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_res, "GBM", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  fit <- train_survivalsvm(est_dd)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "survival - SVM", list_train_vali_Data)
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
    for (alpha in seq(0.1, 0.9, 0.1)) {
      fit <- train_enet(est_dd_rsf, rsf_vars, alpha, seed)
      rs <- calculate_risk_scores(val_dd_list_rsf, function(x) predict_enet(fit, x, rsf_vars))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("RSF + Enet[α=", alpha, "]"), list_train_vali_Data)
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
    rs <- lapply(val_dd_list_rsf, function(x) cbind(x[, 1:2], RS = predict_superpc(spc_r$fit, spc_r$cv_fit, est_dd_rsf, x)))
    tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "RSF + SuperPC", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    fit <- train_survivalsvm(est_dd_rsf)
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
  # PHASE 4: StepCox combinations (51 = 3 directions × 17)
  # ============================================
  message("--- Phase 4: StepCox combinations (51) ---")

  for (direction in c("both", "backward", "forward")) {
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
      for (alpha in seq(0.1, 0.9, 0.1)) {
        fit <- train_enet(est_dd_sc, sc_vars, alpha, seed)
        rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_enet(fit, x, sc_vars))
        tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + Enet[α=", alpha, "]"), list_train_vali_Data)
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
      set.seed(seed)
      Surv <- survival::Surv
      fit <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_dd_sc, ntree = 1000,
                                     nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict(fit, newdata = x)$predicted)
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + RSF"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + SuperPC
      spc_r <- train_superpc(est_dd_sc, seed)
      rs <- lapply(val_dd_list_sc, function(x) cbind(x[, 1:2], RS = predict_superpc(spc_r$fit, spc_r$cv_fit, est_dd_sc, x)))
      tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, paste0(prefix, " + SuperPC"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

      # StepCox + survival-SVM
      fit <- train_survivalsvm(est_dd_sc)
      rs <- calculate_risk_scores(val_dd_list_sc, function(x) predict_survivalsvm(fit, x))
      tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0(prefix, " + survival-SVM"), list_train_vali_Data)
      result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
    }
  }

  # ============================================
  # PHASE 5: CoxBoost combinations (19)
  # ============================================
  message("--- Phase 5: CoxBoost combinations (19) ---")

  # CoxBoost + Enet (9)
  for (alpha in seq(0.1, 0.9, 0.1)) {
    fit <- train_enet(est_dd, pre_var, alpha, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("CoxBoost + Enet[α=", alpha, "]"), list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # CoxBoost + GBM
  gbm_r <- train_gbm(est_dd, seed, cores_for_parallel)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm_r$fit, gbm_r$best, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, gbm_r, "CoxBoost + GBM", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # CoxBoost + Lasso
  fit <- train_lasso(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + Lasso", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # CoxBoost + plsRcox
  fit <- train_plsrcox(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + plsRcox", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # CoxBoost + Ridge
  fit <- train_ridge(est_dd, pre_var, seed)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + Ridge", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # CoxBoost + StepCox (3 directions)
  for (direction in c("both", "backward", "forward")) {
    fit <- train_stepcox(est_dd, direction)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, paste0("CoxBoost + StepCox[", direction, "]"), list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # CoxBoost + SuperPC
  spc_r <- train_superpc(est_dd, seed)
  rs <- lapply(val_dd_list, function(x) cbind(x[, 1:2], RS = predict_superpc(spc_r$fit, spc_r$cv_fit, est_dd, x)))
  tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "CoxBoost + SuperPC", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # CoxBoost + survival-SVM
  fit <- train_survivalsvm(est_dd)
  rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
  tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "CoxBoost + survival-SVM", list_train_vali_Data)
  result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

  # ============================================
  # PHASE 6: Lasso combinations (9)
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
    set.seed(seed)
    Surv <- survival::Surv
    fit <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_dd_lasso, ntree = 1000,
                                   nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict(fit, newdata = x)$predicted)
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
    rs <- lapply(val_dd_list_lasso, function(x) cbind(x[, 1:2], RS = predict_superpc(spc_r$fit, spc_r$cv_fit, est_dd_lasso, x)))
    tmp <- add_model_result(result, ml.res, riskscore, rs, spc_r, "Lasso + SuperPC", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore

    # Lasso + survival-SVM
    fit <- train_survivalsvm(est_dd_lasso)
    rs <- calculate_risk_scores(val_dd_list_lasso, function(x) predict_survivalsvm(fit, x))
    tmp <- add_model_result(result, ml.res, riskscore, rs, fit, "Lasso + survival-SVM", list_train_vali_Data)
    result <- tmp$result; ml.res <- tmp$ml.res; riskscore <- tmp$riskscore
  }

  # ============================================
  # Summary
  # ============================================
  message(paste0("--- Total models: ", length(ml.res), " ---"))

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
  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # Get first algorithm's selected variables
  if (double_ml1 == "RSF") {
    set.seed(seed)
    Surv <- survival::Surv
    fit1 <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000,
                                    nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
    selected_vars <- randomForestSRC::var.select(object = fit1, conservative = "high")$topvars
  } else if (double_ml1 == "StepCox") {
    fit1 <- train_stepcox(est_dd, direction_for_stepcox)
    selected_vars <- get_stepcox_selected_vars(fit1)
  } else if (double_ml1 == "CoxBoost") {
    fit1 <- train_coxboost(est_dd, seed)
    selected_vars <- pre_var  # CoxBoost uses all features
  } else if (double_ml1 == "Lasso") {
    fit1 <- train_lasso(est_dd, pre_var, seed)
    selected_vars <- get_lasso_selected_vars(fit1)
  }

  if (length(selected_vars) <= 1) {
    warning("First algorithm selected < 2 variables")
    return(list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var))
  }

  est_dd2 <- train_data[, c("OS.time", "OS", selected_vars)]
  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
    x[, c("OS.time", "OS", selected_vars)]
  })

  # Run second algorithm
  if (double_ml2 == "CoxBoost") {
    fit2 <- train_coxboost(est_dd2, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_coxboost(fit2, x))
    model_name <- paste0(double_ml1, " + CoxBoost")
  } else if (double_ml2 == "Enet") {
    fit2 <- train_enet(est_dd2, selected_vars, alpha_for_enet, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_enet(fit2, x, selected_vars))
    model_name <- paste0(double_ml1, " + Enet[α=", alpha_for_enet, "]")
  } else if (double_ml2 == "GBM") {
    gbm_res <- train_gbm(est_dd2, seed, cores_for_parallel)
    fit2 <- gbm_res
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_gbm(gbm_res$fit, gbm_res$best, x))
    model_name <- paste0(double_ml1, " + GBM")
  } else if (double_ml2 == "Lasso") {
    fit2 <- train_lasso(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_lasso(fit2, x, selected_vars))
    model_name <- paste0(double_ml1, " + Lasso")
  } else if (double_ml2 == "plsRcox") {
    fit2 <- train_plsrcox(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_plsrcox(fit2, x))
    model_name <- paste0(double_ml1, " + plsRcox")
  } else if (double_ml2 == "Ridge") {
    fit2 <- train_ridge(est_dd2, selected_vars, seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_ridge(fit2, x, selected_vars))
    model_name <- paste0(double_ml1, " + Ridge")
  } else if (double_ml2 == "RSF") {
    set.seed(seed)
    Surv <- survival::Surv
    fit2 <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000,
                                    nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict(fit2, newdata = x)$predicted)
    model_name <- paste0(double_ml1, " + RSF")
  } else if (double_ml2 == "StepCox") {
    fit2 <- train_stepcox(est_dd2, direction_for_stepcox)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_stepcox(fit2, x))
    model_name <- paste0(double_ml1, " + StepCox[", direction_for_stepcox, "]")
  } else if (double_ml2 == "SuperPC") {
    spc_res <- train_superpc(est_dd2, seed)
    fit2 <- spc_res
    rs <- lapply(val_dd_list2, function(x) {
      cbind(x[, 1:2], RS = predict_superpc(spc_res$fit, spc_res$cv_fit, est_dd2, x))
    })
    model_name <- paste0(double_ml1, " + SuperPC")
  } else if (double_ml2 == "survivalsvm") {
    fit2 <- train_survivalsvm(est_dd2)
    rs <- calculate_risk_scores(val_dd_list2, function(x) predict_survivalsvm(fit2, x))
    model_name <- paste0(double_ml1, " + survival-SVM")
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
