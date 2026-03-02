# Parallel Execution for 117 Algorithm Combinations
#
# Linux-only using fork (mclapply)
# Fixed 12 cores with progress display

#' Run all 117 combinations in parallel
#' @inheritParams run_all_algorithms_117
#' @param cores Number of cores (default: 12)
#' @return Combined results
#' @keywords internal
run_all_algorithms_117_parallel <- function(est_dd, train_data, val_dd_list,
                                             list_train_vali_Data, pre_var,
                                             rf_nodesize, seed, cores = 12) {

  result <- data.frame()
  ml.res <- list()
  riskscore <- list()

  # ============================================
  # PHASE 1: Feature selectors (MUST be sequential - shared state)
  # ============================================
  cat("\n=== Phase 1/2: Training feature selectors ===\n")

  set.seed(seed)
  Surv <- survival::Surv

  # RSF
  rsf_fit <- randomForestSRC::rfsrc(
    Surv(OS.time, OS) ~ ., data = est_dd,
    ntree = 1000, nodesize = rf_nodesize,
    splitrule = "logrank", importance = TRUE,
    proximity = TRUE, forest = TRUE, seed = seed
  )
  set.seed(seed)
  rsf_vars <- randomForestSRC::var.select(object = rsf_fit, conservative = "high")$topvars
  cat(sprintf("  RSF: %d vars\n", length(rsf_vars)))

  # StepCox (3 directions)
  stepcox_fits <- list()
  stepcox_vars <- list()
  for (dir in c("both", "backward", "forward")) {
    stepcox_fits[[dir]] <- train_stepcox(est_dd, dir)
    stepcox_vars[[dir]] <- get_stepcox_selected_vars(stepcox_fits[[dir]])
    cat(sprintf("  StepCox[%s]: %d vars\n", dir, length(stepcox_vars[[dir]])))
  }

  # CoxBoost
  coxboost_fit <- train_coxboost(est_dd, seed)
  coxboost_vars <- pre_var

  # Lasso
  lasso_fit <- train_lasso(est_dd, pre_var, seed)
  lasso_vars <- get_lasso_selected_vars(lasso_fit)
  cat(sprintf("  Lasso: %d vars\n\n", length(lasso_vars)))

  # ============================================
  # PHASE 2: Generate ALL parallel tasks
  # ============================================
  cat("=== Phase 2/2: Running 117 combinations in parallel (12 cores) ===\n")

  tasks <- list()

  # ---- Phase 2A: Single models (20) ----
  # RSF
  tasks[["RSF"]] <- function() {
    rs <- calculate_risk_scores(val_dd_list, function(x) predict(rsf_fit, newdata = x)$predicted)
    list(name = "RSF", rs = rs, fit = rsf_fit)
  }

  # Enet (9)
  for (alpha in seq(0.1, 0.9, 0.1)) {
    nm <- paste0("Enet[α=", alpha, "]")
    a <- alpha  # Capture
    tasks[[nm]] <- function() {
      fit <- train_enet(est_dd, pre_var, a, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
      list(name = paste0("Enet[α=", a, "]"), rs = rs, fit = fit)
    }
  }

  # StepCox (3) - use cached fits
  for (dir in c("both", "backward", "forward")) {
    nm <- paste0("StepCox[", dir, "]")
    d <- dir
    tasks[[nm]] <- function() {
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(stepcox_fits[[d]], x))
      list(name = paste0("StepCox[", d, "]"), rs = rs, fit = stepcox_fits[[d]])
    }
  }

  # CoxBoost
  tasks[["CoxBoost"]] <- function() {
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_coxboost(coxboost_fit, x))
    list(name = "CoxBoost", rs = rs, fit = coxboost_fit)
  }

  # Other singles
  tasks[["plsRcox"]] <- function() {
    fit <- train_plsrcox(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
    list(name = "plsRcox", rs = rs, fit = fit)
  }

  tasks[["SuperPC"]] <- function() {
    spc <- train_superpc(est_dd, seed)
    rs <- lapply(val_dd_list, function(x) cbind(x[, 1:2], RS = predict_superpc(spc$fit, spc$cv_fit, est_dd, x)))
    list(name = "SuperPC", rs = rs, fit = spc)
  }

  tasks[["GBM"]] <- function() {
    gbm <- train_gbm(est_dd, seed, cores)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm$fit, gbm$best, x))
    list(name = "GBM", rs = rs, fit = gbm)
  }

  tasks[["survival - SVM"]] <- function() {
    fit <- train_survivalsvm(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
    list(name = "survival - SVM", rs = rs, fit = fit)
  }

  tasks[["Ridge"]] <- function() {
    fit <- train_ridge(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
    list(name = "Ridge", rs = rs, fit = fit)
  }

  # Lasso
  tasks[["Lasso"]] <- function() {
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(lasso_fit, x, pre_var))
    list(name = "Lasso", rs = rs, fit = lasso_fit)
  }

  # ---- Phase 2B: RSF combinations (19) ----
  if (length(rsf_vars) > 1) {
    est_rsf <- train_data[, c("OS.time", "OS", rsf_vars)]
    val_rsf <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", rsf_vars)])

    # RSF + CoxBoost
    tasks[["RSF + CoxBoost"]] <- function() {
      fit <- train_coxboost(est_rsf, seed)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_coxboost(fit, x))
      list(name = "RSF + CoxBoost", rs = rs, fit = fit)
    }

    # RSF + Enet (9)
    for (alpha in seq(0.1, 0.9, 0.1)) {
      nm <- paste0("RSF + Enet[α=", alpha, "]")
      a <- alpha
      tasks[[nm]] <- function() {
        fit <- train_enet(est_rsf, rsf_vars, a, seed)
        rs <- calculate_risk_scores(val_rsf, function(x) predict_enet(fit, x, rsf_vars))
        list(name = paste0("RSF + Enet[α=", a, "]"), rs = rs, fit = fit)
      }
    }

    # RSF + others
    tasks[["RSF + GBM"]] <- function() {
      gbm <- train_gbm(est_rsf, seed, cores)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_gbm(gbm$fit, gbm$best, x))
      list(name = "RSF + GBM", rs = rs, fit = gbm)
    }

    tasks[["RSF + Lasso"]] <- function() {
      fit <- train_lasso(est_rsf, rsf_vars, seed)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_lasso(fit, x, rsf_vars))
      list(name = "RSF + Lasso", rs = rs, fit = fit)
    }

    tasks[["RSF + plsRcox"]] <- function() {
      fit <- train_plsrcox(est_rsf, rsf_vars, seed)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_plsrcox(fit, x))
      list(name = "RSF + plsRcox", rs = rs, fit = fit)
    }

    tasks[["RSF + Ridge"]] <- function() {
      fit <- train_ridge(est_rsf, rsf_vars, seed)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_ridge(fit, x, rsf_vars))
      list(name = "RSF + Ridge", rs = rs, fit = fit)
    }

    tasks[["RSF + SuperPC"]] <- function() {
      spc <- train_superpc(est_rsf, seed)
      rs <- lapply(val_rsf, function(x) cbind(x[, 1:2], RS = predict_superpc(spc$fit, spc$cv_fit, est_rsf, x)))
      list(name = "RSF + SuperPC", rs = rs, fit = spc)
    }

    tasks[["RSF + survival-SVM"]] <- function() {
      fit <- train_survivalsvm(est_rsf, seed)
      rs <- calculate_risk_scores(val_rsf, function(x) predict_survivalsvm(fit, x))
      list(name = "RSF + survival-SVM", rs = rs, fit = fit)
    }

    # RSF + StepCox (3)
    for (dir in c("both", "backward", "forward")) {
      nm <- paste0("RSF + StepCox[", dir, "]")
      d <- dir
      tasks[[nm]] <- function() {
        fit <- train_stepcox(est_rsf, d)
        rs <- calculate_risk_scores(val_rsf, function(x) predict_stepcox(fit, x))
        list(name = paste0("RSF + StepCox[", d, "]"), rs = rs, fit = fit)
      }
    }
  }

  # ---- Phase 2C: StepCox combinations (51 = 3 × 17) ----
  for (dir in c("both", "backward", "forward")) {
    sc_vars <- stepcox_vars[[dir]]
    if (length(sc_vars) > 1) {
      est_sc <- train_data[, c("OS.time", "OS", sc_vars)]
      val_sc <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", sc_vars)])
      prefix <- paste0("StepCox[", dir, "]")
      d <- dir

      tasks[[paste0(prefix, " + CoxBoost")]] <- function() {
        fit <- train_coxboost(est_sc, seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict_coxboost(fit, x))
        list(name = paste0("StepCox[", d, "] + CoxBoost"), rs = rs, fit = fit)
      }

      for (alpha in seq(0.1, 0.9, 0.1)) {
        nm <- paste0(prefix, " + Enet[α=", alpha, "]")
        a <- alpha
        tasks[[nm]] <- function() {
          fit <- train_enet(est_sc, sc_vars, a, seed)
          rs <- calculate_risk_scores(val_sc, function(x) predict_enet(fit, x, sc_vars))
          list(name = paste0("StepCox[", d, "] + Enet[α=", a, "]"), rs = rs, fit = fit)
        }
      }

      tasks[[paste0(prefix, " + GBM")]] <- function() {
        gbm <- train_gbm(est_sc, seed, cores)
        rs <- calculate_risk_scores(val_sc, function(x) predict_gbm(gbm$fit, gbm$best, x))
        list(name = paste0("StepCox[", d, "] + GBM"), rs = rs, fit = gbm)
      }

      tasks[[paste0(prefix, " + Lasso")]] <- function() {
        fit <- train_lasso(est_sc, sc_vars, seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict_lasso(fit, x, sc_vars))
        list(name = paste0("StepCox[", d, "] + Lasso"), rs = rs, fit = fit)
      }

      tasks[[paste0(prefix, " + plsRcox")]] <- function() {
        fit <- train_plsrcox(est_sc, sc_vars, seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict_plsrcox(fit, x))
        list(name = paste0("StepCox[", d, "] + plsRcox"), rs = rs, fit = fit)
      }

      tasks[[paste0(prefix, " + Ridge")]] <- function() {
        fit <- train_ridge(est_sc, sc_vars, seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict_ridge(fit, x, sc_vars))
        list(name = paste0("StepCox[", d, "] + Ridge"), rs = rs, fit = fit)
      }

      tasks[[paste0(prefix, " + RSF")]] <- function() {
        fit <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_sc, ntree = 1000,
                                       nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict(fit, newdata = x)$predicted)
        list(name = paste0("StepCox[", d, "] + RSF"), rs = rs, fit = fit)
      }

      tasks[[paste0(prefix, " + SuperPC")]] <- function() {
        spc <- train_superpc(est_sc, seed)
        rs <- lapply(val_sc, function(x) cbind(x[, 1:2], RS = predict_superpc(spc$fit, spc$cv_fit, est_sc, x)))
        list(name = paste0("StepCox[", d, "] + SuperPC"), rs = rs, fit = spc)
      }

      tasks[[paste0(prefix, " + survival-SVM")]] <- function() {
        fit <- train_survivalsvm(est_sc, seed)
        rs <- calculate_risk_scores(val_sc, function(x) predict_survivalsvm(fit, x))
        list(name = paste0("StepCox[", d, "] + survival-SVM"), rs = rs, fit = fit)
      }
    }
  }

  # ---- Phase 2D: CoxBoost combinations (19) ----
  # CoxBoost + Enet (9)
  for (alpha in seq(0.1, 0.9, 0.1)) {
    nm <- paste0("CoxBoost + Enet[α=", alpha, "]")
    a <- alpha
    tasks[[nm]] <- function() {
      fit <- train_enet(est_dd, pre_var, a, seed)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_enet(fit, x, pre_var))
      list(name = paste0("CoxBoost + Enet[α=", a, "]"), rs = rs, fit = fit)
    }
  }

  tasks[["CoxBoost + GBM"]] <- function() {
    gbm <- train_gbm(est_dd, seed, cores)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_gbm(gbm$fit, gbm$best, x))
    list(name = "CoxBoost + GBM", rs = rs, fit = gbm)
  }

  tasks[["CoxBoost + Lasso"]] <- function() {
    fit <- train_lasso(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_lasso(fit, x, pre_var))
    list(name = "CoxBoost + Lasso", rs = rs, fit = fit)
  }

  tasks[["CoxBoost + plsRcox"]] <- function() {
    fit <- train_plsrcox(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_plsrcox(fit, x))
    list(name = "CoxBoost + plsRcox", rs = rs, fit = fit)
  }

  tasks[["CoxBoost + Ridge"]] <- function() {
    fit <- train_ridge(est_dd, pre_var, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_ridge(fit, x, pre_var))
    list(name = "CoxBoost + Ridge", rs = rs, fit = fit)
  }

  # CoxBoost + StepCox (3)
  for (dir in c("both", "backward", "forward")) {
    nm <- paste0("CoxBoost + StepCox[", dir, "]")
    d <- dir
    tasks[[nm]] <- function() {
      fit <- train_stepcox(est_dd, d)
      rs <- calculate_risk_scores(val_dd_list, function(x) predict_stepcox(fit, x))
      list(name = paste0("CoxBoost + StepCox[", d, "]"), rs = rs, fit = fit)
    }
  }

  tasks[["CoxBoost + SuperPC"]] <- function() {
    spc <- train_superpc(est_dd, seed)
    rs <- lapply(val_dd_list, function(x) cbind(x[, 1:2], RS = predict_superpc(spc$fit, spc$cv_fit, est_dd, x)))
    list(name = "CoxBoost + SuperPC", rs = rs, fit = spc)
  }

  tasks[["CoxBoost + survival-SVM"]] <- function() {
    fit <- train_survivalsvm(est_dd, seed)
    rs <- calculate_risk_scores(val_dd_list, function(x) predict_survivalsvm(fit, x))
    list(name = "CoxBoost + survival-SVM", rs = rs, fit = fit)
  }

  # ---- Phase 2E: Lasso combinations (9) ----
  if (length(lasso_vars) > 1) {
    est_lasso <- train_data[, c("OS.time", "OS", lasso_vars)]
    val_lasso <- lapply(list_train_vali_Data, function(x) x[, c("OS.time", "OS", lasso_vars)])

    tasks[["Lasso + CoxBoost"]] <- function() {
      fit <- train_coxboost(est_lasso, seed)
      rs <- calculate_risk_scores(val_lasso, function(x) predict_coxboost(fit, x))
      list(name = "Lasso + CoxBoost", rs = rs, fit = fit)
    }

    tasks[["Lasso + GBM"]] <- function() {
      gbm <- train_gbm(est_lasso, seed, cores)
      rs <- calculate_risk_scores(val_lasso, function(x) predict_gbm(gbm$fit, gbm$best, x))
      list(name = "Lasso + GBM", rs = rs, fit = gbm)
    }

    tasks[["Lasso + plsRcox"]] <- function() {
      fit <- train_plsrcox(est_lasso, lasso_vars, seed)
      rs <- calculate_risk_scores(val_lasso, function(x) predict_plsrcox(fit, x))
      list(name = "Lasso + plsRcox", rs = rs, fit = fit)
    }

    tasks[["Lasso + RSF"]] <- function() {
      fit <- randomForestSRC::rfsrc(Surv(OS.time, OS) ~ ., data = est_lasso, ntree = 1000,
                                     nodesize = rf_nodesize, splitrule = "logrank", seed = seed)
      rs <- calculate_risk_scores(val_lasso, function(x) predict(fit, newdata = x)$predicted)
      list(name = "Lasso + RSF", rs = rs, fit = fit)
    }

    # Lasso + StepCox (3)
    for (dir in c("both", "backward", "forward")) {
      nm <- paste0("Lasso + StepCox[", dir, "]")
      d <- dir
      tasks[[nm]] <- function() {
        fit <- train_stepcox(est_lasso, d)
        rs <- calculate_risk_scores(val_lasso, function(x) predict_stepcox(fit, x))
        list(name = paste0("Lasso + StepCox[", d, "]"), rs = rs, fit = fit)
      }
    }

    tasks[["Lasso + SuperPC"]] <- function() {
      spc <- train_superpc(est_lasso, seed)
      rs <- lapply(val_lasso, function(x) cbind(x[, 1:2], RS = predict_superpc(spc$fit, spc$cv_fit, est_lasso, x)))
      list(name = "Lasso + SuperPC", rs = rs, fit = spc)
    }

    tasks[["Lasso + survival-SVM"]] <- function() {
      fit <- train_survivalsvm(est_lasso, seed)
      rs <- calculate_risk_scores(val_lasso, function(x) predict_survivalsvm(fit, x))
      list(name = "Lasso + survival-SVM", rs = rs, fit = fit)
    }
  }

  # ============================================
  # Execute tasks in parallel with progress
  # ============================================
  n_tasks <- length(tasks)
  cat(sprintf("  Total tasks: %d\n", n_tasks))

  pb <- txtProgressBar(min = 0, max = n_tasks, style = 3, width = 60)

  results <- parallel::mclapply(
    names(tasks),
    function(nm) {
      res <- tryCatch(
        tasks[[nm]](),
        error = function(e) list(name = nm, error = conditionMessage(e))
      )
      # Update progress (note: may not display in parallel)
      res
    },
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = TRUE
  )

  # Show final progress
  setTxtProgressBar(pb, n_tasks)
  close(pb)

  names(results) <- names(tasks)

  # ============================================
  # Combine results
  # ============================================
  cat("\n  Combining results...\n")

  for (nm in names(results)) {
    r <- results[[nm]]
    if (!is.null(r$error)) {
      cat(sprintf("  [WARN] %s: %s\n", nm, r$error))
      next
    }

    rs <- return_id_to_rs(r$rs, list_train_vali_Data)
    cc <- calculate_cindex_result(rs, r$name)
    result <- rbind(result, cc)
    ml.res[[r$name]] <- r$fit
    riskscore[[r$name]] <- rs
  }

  cat(sprintf("\n=== Complete: %d models built ===\n", length(ml.res)))

  return(list(
    Cindex.res = result,
    ml.res = ml.res,
    riskscore = riskscore,
    Sig.genes = pre_var
  ))
}
