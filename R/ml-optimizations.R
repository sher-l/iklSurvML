# Performance Optimizations for Mime Package
#
# This file contains optimized versions of key functions with:
# 1. Caching for repeated computations
# 2. Parallel execution of independent algorithms
# 3. Progress tracking

# ---- Caching System ----

#' Create a simple cache environment
#' @keywords internal
ml_cache <- new.env(parent = emptyenv())

#' Get cached result or compute
#' @param key Cache key
#' @param compute_fn Function to compute if not cached
#' @return Cached or computed result
#' @keywords internal
cached_compute <- function(key, compute_fn) {
  if (exists(key, envir = ml_cache)) {
    message(paste0("Using cached result for: ", key))
    return(get(key, envir = ml_cache))
  }
  result <- compute_fn()
  assign(key, result, envir = ml_cache)
  return(result)
}

#' Clear the ML cache
#' @export
clear_ml_cache <- function() {
  rm(list = ls(envir = ml_cache), envir = ml_cache)
  message("ML cache cleared")
}

# ---- Parallel Algorithm Execution ----

#' Run all algorithms in parallel using future
#' @param est_dd Training data
#' @param train_data Full training data
#' @param val_dd_list Validation data list
#' @param list_train_vali_Data Original data list
#' @param pre_var Predictor variables
#' @param rf_nodesize RSF node size
#' @param seed Random seed
#' @param cores_for_parallel Number of cores
#' @return Combined results
#' @keywords internal
run_all_algorithms_parallel <- function(est_dd,
                                        train_data,
                                        val_dd_list,
                                        list_train_vali_Data,
                                        pre_var,
                                        rf_nodesize = 5,
                                        seed = 5201314,
                                        cores_for_parallel = 6) {

  # Set up parallel backend
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession, workers = cores_for_parallel)
    on.exit(future::plan(future::sequential))
    use_parallel <- TRUE
  } else {
    warning("future/future.apply not available, running sequentially")
    use_parallel <- FALSE
  }

  # Define algorithm tasks (these can run in parallel)
  algorithm_tasks <- list(
    # Single algorithms
    RSF = function() {
      set.seed(seed)
      fit <- train_rsf(est_dd, rf_nodesize, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_rsf(fit, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "RSF")
    },

    CoxBoost = function() {
      set.seed(seed)
      fit <- train_coxboost(est_dd, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(fit, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "CoxBoost")
    },

    plsRcox = function() {
      fit <- train_plsrcox(est_dd, pre_var, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "plsRcox")
    },

    Ridge = function() {
      fit <- train_ridge(est_dd, pre_var, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "Ridge")
    },

    Lasso = function() {
      fit <- train_lasso(est_dd, pre_var, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "Lasso")
    },

    survivalsvm = function() {
      fit <- train_survivalsvm(est_dd, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = fit, rs = rs, name = "survival-SVM")
    },

    GBM = function() {
      gbm_result <- train_gbm(est_dd, seed, cores_for_parallel)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm_result$fit, gbm_result$best, x))
      rs <- return_id_to_rs(rs, list_train_vali_Data)
      list(model = gbm_result, rs = rs, name = "GBM")
    }
  )

  # Run algorithms in parallel
  if (use_parallel) {
    results <- future.apply::future_lapply(algorithm_tasks, function(f) f())
  } else {
    results <- lapply(algorithm_tasks, function(f) f())
  }

  # Combine results
  result_df <- data.frame()
  ml.res <- list()
  riskscore <- list()

  for (name in names(results)) {
    r <- results[[name]]
    cc <- calculate_cindex_result(r$rs, r$name)
    result_df <- rbind(result_df, cc)
    ml.res[[r$name]] <- r$model
    riskscore[[r$name]] <- r$rs
  }

  return(list(
    Cindex.res = result_df,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}

# ---- Optimized RSF with Caching ----

#' Train RSF with caching for feature selection
#' @inheritParams train_rsf
#' @return List with fit and selected variables
#' @keywords internal
train_rsf_cached <- function(est_dd, rf_nodesize = 5, seed = 5201314) {
  cache_key <- paste0("rsf_", digest::digest(list(
    nrow(est_dd), ncol(est_dd), rf_nodesize, seed
  )))

  cached_compute(cache_key, function() {
    fit <- train_rsf(est_dd, rf_nodesize, seed)
    rid <- get_rsf_selected_vars(fit, seed)
    list(fit = fit, selected_vars = rid)
  })
}

# ---- Batch Processing ----

#' Process algorithms in batches to manage memory
#' @param tasks List of algorithm tasks
#' @param batch_size Number of tasks per batch
#' @return Combined results
#' @keywords internal
run_in_batches <- function(tasks, batch_size = 4) {
  n_tasks <- length(tasks)
  n_batches <- ceiling(n_tasks / batch_size)

  all_results <- list()

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_tasks)
    batch_tasks <- tasks[start_idx:end_idx]

    message(paste0("Processing batch ", i, "/", n_batches, " (tasks ", start_idx, "-", end_idx, ")"))

    batch_results <- lapply(batch_tasks, function(f) f())
    all_results <- c(all_results, batch_results)

    # Force garbage collection between batches
    gc()
  }

  return(all_results)
}

# ---- Progress Tracking ----

#' Create progress bar for long operations
#' @param total Total number of iterations
#' @return Progress function
#' @keywords internal
create_progress_tracker <- function(total) {
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = "  [:bar] :percent eta: :eta",
      total = total,
      clear = FALSE,
      width = 60
    )
    return(function() pb$tick())
  } else {
    return(function() {
      cat(".")
    })
  }
}
