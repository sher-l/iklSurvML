make_survival_smoke_data <- function(n = 30, p = 5, prefix = "S") {
  idx <- seq_len(n * p)
  x <- matrix(sin(idx * 0.37) + cos(idx * 0.13), nrow = n, ncol = p)
  colnames(x) <- paste0("G", seq_len(p))
  data.frame(
    ID = paste0(prefix, seq_len(n)),
    OS.time = seq_len(n) + 1,
    OS = rep(c(1, 0, 1), length.out = n),
    x,
    check.names = FALSE
  )
}

make_mock_all_mode_result <- function(n = 117, prefix = "Mock") {
  model_names <- paste0(prefix, sprintf("%03d", seq_len(n)))
  fits <- stats::setNames(
    lapply(seq_len(n), function(i) structure(list(index = i), class = "mock_survival_fit")),
    model_names
  )
  riskscores <- stats::setNames(
    lapply(seq_len(n), function(i) {
      data.frame(
        ID = "Train1",
        OS.time = 1,
        OS = 1,
        RS = i,
        stringsAsFactors = FALSE
      )
    }),
    model_names
  )

  list(
    Cindex.res = data.frame(
      ID = "Train",
      Cindex = rep(0.5, n),
      Model = model_names,
      stringsAsFactors = FALSE
    ),
    ml.res = fits,
    riskscore = riskscores,
    Sig.genes = paste0("G", seq_len(4))
  )
}

test_that("calculate_cindex preserves risk-score direction without refitting", {
  data <- data.frame(
    OS.time = seq_len(10),
    OS = rep(1, 10)
  )

  high_risk_dies_earlier <- -data$OS.time
  high_risk_dies_later <- data$OS.time

  expect_equal(iklSurvML:::calculate_cindex(high_risk_dies_earlier, data), 1)
  expect_equal(iklSurvML:::calculate_cindex(high_risk_dies_later, data), 0)
})

test_that("ML.Dev.Prog.Sig.Fast uses the documented default seed", {
  skip_if_not_installed("glmnet")

  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "T")
  validation_data <- make_survival_smoke_data(n = 25, p = 6, prefix = "V")
  candidate_genes <- paste0("G", seq_len(6))

  result <- NULL
  expect_no_error({
    result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data, Val = validation_data),
      candidate_genes = candidate_genes,
      unicox.filter.for.candi = FALSE,
      mode = "single",
      single_ml = "Lasso",
      nodesize = 5,
      use_parallel = FALSE
    )
  })

  expect_equal(names(result$ml.res), "Lasso")
})

test_that("ML.Dev.Prog.Sig.Fast defaults to the fork-safe sequential all-mode", {
  expect_false(formals(iklSurvML::ML.Dev.Prog.Sig.Fast)$use_parallel)
})

test_that("ML.Dev.Prog.Sig.Fast all-mode defaults to the complete 117-grid result", {
  train_data <- make_survival_smoke_data(n = 12, p = 4, prefix = "F")

  testthat::local_mocked_bindings(
    run_all_algorithms_128 = function(...) {
      expect_equal(list(...)$model_grid, "117")
      make_mock_all_mode_result(117L, "FastMock")
    },
    .package = "iklSurvML"
  )

  result <- NULL
  expect_no_error({
    result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      unicox.filter.for.candi = FALSE,
      mode = "all",
      use_parallel = FALSE
    )
  })

  expect_length(result$ml.res, 117)
  expect_length(result$Model.info, 117)
  expect_equal(names(result$ml.res)[[1]], "FastMock001")
})

test_that("ML.Dev.Prog.Sig all-mode uses the fixed 117-grid result", {
  train_data <- make_survival_smoke_data(n = 12, p = 4, prefix = "G")

  testthat::local_mocked_bindings(
    run_all_algorithms = function(...) {
      expect_equal(list(...)$model_grid, "117")
      make_mock_all_mode_result(117L, "Grid117Mock")
    },
    .package = "iklSurvML"
  )

  result <- NULL
  expect_no_error({
    result <- iklSurvML::ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      unicox.filter.for.candi = FALSE,
      mode = "all"
    )
  })

  expect_length(result$ml.res, 117)
  expect_length(result$Model.info, 117)
  expect_equal(names(result$ml.res)[[1]], "Grid117Mock001")
})

test_that("all-mode entry points reject partial grids unless explicitly allowed", {
  train_data <- make_survival_smoke_data(n = 12, p = 4, prefix = "P")

  testthat::local_mocked_bindings(
    run_all_algorithms = function(...) make_mock_all_mode_result(100L, "PartialMock"),
    .package = "iklSurvML"
  )

  expect_error(
    iklSurvML::ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      unicox.filter.for.candi = FALSE,
      mode = "all"
    ),
    "allow_partial = TRUE"
  )

  result <- NULL
  expect_no_error({
    result <- iklSurvML::ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      unicox.filter.for.candi = FALSE,
      mode = "all",
      allow_partial = TRUE
    )
  })
  expect_length(result$ml.res, 100)
})

test_that("feature-name normalization handles underscores during unicox filtering", {
  skip_if_not_installed("glmnet")

  train_data <- make_survival_smoke_data(n = 40, p = 4, prefix = "U")
  colnames(train_data)[4:7] <- c("Gene_A", "Gene_B", "Gene_C", "Gene_D")

  result <- NULL
  expect_no_error({
    result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = c("Gene_A", "Gene_B", "Gene_C", "Gene_D"),
      unicox.filter.for.candi = TRUE,
      unicox_p_cutoff = 1,
      mode = "single",
      single_ml = "Lasso",
      use_parallel = FALSE
    )
  })

  expect_equal(names(result$ml.res), "Lasso")
  expect_true(all(grepl("^Gene\\.", result$Sig.genes)))
})

test_that("training data receives the same NA imputation as validation data", {
  skip_if_not_installed("glmnet")

  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "N")
  train_data$G2[5] <- NA_real_

  result <- NULL
  expect_warning(
    expect_warning({
      result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
        train_data = train_data,
        list_train_vali_Data = list(Train = train_data),
        candidate_genes = paste0("G", seq_len(6)),
        unicox.filter.for.candi = FALSE,
        mode = "single",
        single_ml = "Lasso",
        use_parallel = FALSE
      )
    },
    "Dataset 'Train'.*imputed with column means"
  ),
    "Training data.*imputed with column means"
  )

  expect_equal(names(result$ml.res), "Lasso")
})

test_that("survival preprocessing rejects duplicated sample IDs", {
  train_data <- make_survival_smoke_data(n = 10, p = 3, prefix = "DUP")
  train_data$ID[2] <- train_data$ID[1]
  common_feature <- c("ID", "OS.time", "OS", paste0("G", seq_len(3)))

  expect_error(
    iklSurvML:::preprocess_train_data(train_data, common_feature),
    "duplicated sample IDs"
  )
})

test_that("survival preprocessing recipe imputes validation data with training means", {
  train_data <- make_survival_smoke_data(n = 10, p = 3, prefix = "R")
  val_data <- make_survival_smoke_data(n = 8, p = 3, prefix = "V")
  train_data$G1 <- seq_len(nrow(train_data))
  val_data$G1[2] <- NA_real_
  common_feature <- c("ID", "OS.time", "OS", paste0("G", seq_len(3)))

  prepped <- iklSurvML:::preprocess_train_data(
    train_data,
    common_feature,
    return_recipe = TRUE
  )
  val_prepped <- NULL
  expect_warning(
    val_prepped <- iklSurvML:::preprocess_data_list(
      list(Val = val_data),
      common_feature,
      recipe = prepped$recipe
    ),
    "training preprocessing recipe"
  )

  expect_equal(val_prepped$Val$G1[2], mean(train_data$G1))
  expect_equal(prepped$recipe$impute_values[["G1"]], mean(train_data$G1))
})

test_that("survival preprocessing rejects ambiguous event coding", {
  train_data <- make_survival_smoke_data(n = 10, p = 3, prefix = "E")
  train_data$OS <- ifelse(train_data$OS == 1, 2, 1)
  common_feature <- c("ID", "OS.time", "OS", paste0("G", seq_len(3)))

  expect_error(
    iklSurvML:::preprocess_train_data(train_data, common_feature),
    "OS must be coded as 0 = censored/alive and 1 = event/dead"
  )
})

test_that("numeric preprocessing preserves factor labels and rejects nonnumeric labels", {
  train_data <- make_survival_smoke_data(n = 10, p = 3, prefix = "F")
  train_data$G1 <- factor(as.character(seq_len(nrow(train_data))))
  common_feature <- c("ID", "OS.time", "OS", paste0("G", seq_len(3)))

  prepped <- iklSurvML:::preprocess_train_data(train_data, common_feature)
  expect_equal(prepped$G1, as.numeric(seq_len(nrow(train_data))))

  train_data$G1 <- factor(rep(c("low", "high"), length.out = nrow(train_data)))
  expect_error(
    iklSurvML:::preprocess_train_data(train_data, common_feature),
    "contains non-numeric values"
  )
})

test_that("strict feature alignment does not let validation cohorts shrink training features", {
  train_data <- make_survival_smoke_data(n = 10, p = 3, prefix = "T")
  validation_data <- train_data[, setdiff(colnames(train_data), "G3")]

  expect_error(
    iklSurvML:::resolve_survival_common_features(
      train_data = train_data,
      list_train_vali_Data = list(Val = validation_data),
      candidate_genes = paste0("G", seq_len(3)),
      feature_alignment = "strict"
    ),
    "feature_alignment='strict'"
  )

  expect_warning(
    common <- iklSurvML:::resolve_survival_common_features(
      train_data = train_data,
      list_train_vali_Data = list(Val = validation_data),
      candidate_genes = paste0("G", seq_len(3)),
      feature_alignment = "intersection"
    ),
    "dropped candidate features"
  )
  expect_false("G3" %in% common)
})

test_that("survival AUC helper does not flip marker direction using validation outcomes", {
  skip_if_not_installed("survivalROC")

  risk_table <- data.frame(
    OS.time = seq(100, 1000, by = 100),
    OS = rep(1, 10),
    RS = seq_len(10)
  )
  direct <- survivalROC::survivalROC(
    Stime = risk_table$OS.time,
    status = risk_table$OS,
    marker = risk_table$RS,
    predict.time = 365,
    method = "KM"
  )

  observed <- iklSurvML:::calculate_survival_roc_from_risk(
    risk_table,
    AUC_time = 1,
    auc_cal_method = "KM"
  )

  expect_equal(unique(observed$AUC), direct$AUC)
  expect_true(all(observed$marker_direction == "higher_is_worse"))
})

test_that("invalid survival ML parameters fail fast with actionable messages", {
  train_data <- make_survival_smoke_data(n = 30, p = 4, prefix = "P")

  expect_error(
    iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      unicox.filter.for.candi = FALSE,
      mode = "single",
      single_ml = "TypoModel",
      use_parallel = FALSE
    ),
    "single_ml"
  )
})

test_that("survival-SVM model names are canonicalized while accepting aliases", {
  expect_equal(iklSurvML:::normalize_survival_ml_name("survival - SVM"), "survivalsvm")
  expect_equal(iklSurvML:::display_survival_ml_name("survivalsvm"), "survival-SVM")
})

test_that("survival-SVM predictions are oriented as higher risk is worse", {
  skip_if_not_installed("survivalsvm")

  set.seed(777)
  n <- 45
  g1 <- stats::rnorm(n)
  g2 <- stats::rnorm(n)
  dat <- data.frame(
    OS.time = pmax(1, round(1000 - 150 * g1 + stats::rnorm(n, sd = 20))),
    OS = as.integer(stats::runif(n) < stats::plogis(g1)),
    G1 = g1,
    G2 = g2,
    check.names = FALSE
  )

  fit <- iklSurvML:::train_survivalsvm(dat, seed = 123)
  rs <- iklSurvML:::predict_survivalsvm(fit, dat)

  expect_lt(stats::cor(rs, dat$OS.time), 0)
  expect_gt(stats::cor(rs, dat$G1), 0)
})

test_that("StepCox forward selection starts from the null model", {
  set.seed(778)
  n <- 80
  g1 <- stats::rnorm(n)
  dat <- data.frame(
    OS.time = pmax(1, round(900 - 180 * g1 + stats::rnorm(n, sd = 25))),
    OS = rep(1, n),
    G1 = g1,
    G2 = stats::rnorm(n),
    G3 = stats::rnorm(n),
    check.names = FALSE
  )

  full <- survival::coxph(survival::Surv(OS.time, OS) ~ ., dat)
  forward <- iklSurvML:::train_stepcox(dat, "forward")

  expect_true("G1" %in% names(stats::coef(forward)))
  expect_lt(length(stats::coef(forward)), length(stats::coef(full)))
})

test_that("CoxBoost combinations use non-zero coefficients for feature selection", {
  fit <- list(coefficients = c(G1 = 0, G2 = 0.25, G3 = -0.5, G4 = 0))

  expect_equal(iklSurvML:::get_coxboost_selected_vars(fit), c("G2", "G3"))
})

test_that("double-mode rejects self-combinations", {
  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "S")

  expect_error(
    iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(6)),
      unicox.filter.for.candi = FALSE,
      mode = "double",
      double_ml1 = "Lasso",
      double_ml2 = "Lasso",
      use_parallel = FALSE
    ),
    "Self-combinations are not supported"
  )

  expect_error(
    iklSurvML::ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(6)),
      unicox.filter.for.candi = FALSE,
      mode = "double",
      double_ml1 = "CoxBoost",
      double_ml2 = "CoxBoost",
      cores_for_parallel = 1
    ),
    "Self-combinations are not supported"
  )

  expect_error(
    iklSurvML:::validate_survival_ml_params(
      mode = "double",
      double_ml1 = "StepCox",
      double_ml2 = "StepCox",
      alpha_for_enet = 0.1,
      direction_for_stepcox = "both"
    ),
    "Self-combinations are not supported"
  )
})

test_that("double-mode downstream risk and AUC functions use the materialized non-self combination model", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survivalROC")

  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "D")
  train_data$OS.time <- seq(100, by = 100, length.out = nrow(train_data))
  candidate_genes <- paste0("G", seq_len(6))
  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  processed_train <- iklSurvML:::preprocess_train_data(train_data, common_feature)
  est_dd <- as.data.frame(processed_train)[, common_feature[-1]]
  fit <- iklSurvML:::train_lasso(est_dd, candidate_genes, seed = 5201314)
  result <- list(
    ml.res = stats::setNames(list(fit), "RSF + Lasso"),
    Sig.genes = candidate_genes
  )

  rs <- NULL
  expect_no_error({
    rs <- iklSurvML::cal_RS_ml_res(
      res.by.ML.Dev.Prog.Sig = result,
      train_data = train_data,
      inputmatrix.list = list(Train = train_data),
      mode = "double",
      double_ml1 = "RSF",
      double_ml2 = "Lasso"
    )
  })
  expect_equal(names(rs), "RSF + Lasso")

  auc <- NULL
  expect_no_error({
    auc <- iklSurvML::cal_AUC_ml_res(
      res.by.ML.Dev.Prog.Sig = result,
      train_data = train_data,
      inputmatrix.list = list(Train = train_data),
      mode = "double",
      double_ml1 = "RSF",
      double_ml2 = "Lasso",
      AUC_time = 1,
      auc_cal_method = "KM"
    )
  })
  expect_equal(names(auc), "RSF + Lasso")
})

test_that("parallel task errors are promoted instead of silently dropping models", {
  results <- list(
    ok = list(name = "ok", rs = list(), fit = structure(list(), class = "mock")),
    bad = list(name = "bad", error = "boom")
  )

  expect_error(
    iklSurvML:::assert_no_parallel_task_errors(results),
    "bad: boom"
  )
})

test_that("parallel all-mode keeps GBM work out of fork workers", {
  split <- iklSurvML:::split_gbm_tasks_for_fork_safety(c(
    "GBM",
    "RSF + GBM",
    "StepCox[both] + Lasso",
    "Lasso + survival-SVM"
  ))

  expect_equal(split$sequential, c("GBM", "RSF + GBM"))
  expect_equal(split$parallel, c("StepCox[both] + Lasso", "Lasso + survival-SVM"))
})

test_that("GBM training avoids the native cross-validation path that can segfault", {
  body_text <- paste(deparse(iklSurvML:::train_gbm), collapse = "\n")

  expect_true(grepl("cv.folds = 0", body_text, fixed = TRUE))
  expect_false(grepl("fit$cv.error", body_text, fixed = TRUE))
})

test_that("plsRcox component selection retries unstable cross-validation folds", {
  prep <- list(x = data.frame(G1 = seq_len(12), G2 = rev(seq_len(12))))
  calls <- 0L
  fake_cv <- function(...) {
    calls <<- calls + 1L
    if (calls == 1L) {
      stop("unstable fold")
    }
    list(NULL, NULL, NULL, NULL, 2L)
  }

  selected_nt <- iklSurvML:::select_plsrcox_component_count(
    prep = prep,
    time = seq_len(12),
    status = rep(c(1, 0), 6),
    seed = 5201314,
    max_nt = 5,
    cv_fun = fake_cv
  )

  expect_equal(selected_nt, 2L)
  expect_equal(calls, 2L)
})

test_that("SuperPC prediction params use the best CV component and threshold", {
  cv_fit <- list(
    thresholds = c(0.1, 0.2, 0.3),
    scor = matrix(
      c(1, 2, 3,
        4, 9, 6),
      nrow = 2,
      byrow = TRUE
    )
  )

  params <- iklSurvML:::select_superpc_prediction_params(cv_fit)

  expect_equal(params$n.components, 2L)
  expect_equal(params$threshold, 0.2)
})

test_that("SuperPC model wrappers preserve their training feature subset", {
  training <- data.frame(
    OS.time = c(10, 20, 30),
    OS = c(1, 0, 1),
    G1 = c(1, 2, 3),
    G2 = c(3, 2, 1),
    check.names = FALSE
  )

  model <- iklSurvML:::make_superpc_model(
    fit = structure(list(), class = "mock_superpc_fit"),
    cv_fit = list(thresholds = 1, scor = matrix(1, nrow = 1)),
    train_data = training
  )
  parts <- iklSurvML:::extract_superpc_model(model)

  expect_equal(parts$features, c("G1", "G2"))
  expect_equal(colnames(parts$train_data), c("OS.time", "OS", "G1", "G2"))
})

test_that("SuperPC predictions return exactly one risk score per sample", {
  skip_if_not_installed("superpc")

  train_data <- make_survival_smoke_data(n = 45, p = 6, prefix = "SP")
  est_dd <- train_data[, -1]
  fit <- suppressWarnings(iklSurvML:::train_superpc(est_dd, seed = 5201314))
  rs <- suppressWarnings(iklSurvML:::predict_superpc_model(fit, est_dd, est_dd))

  expect_equal(length(rs), nrow(est_dd))
})

test_that("cal_RS_ml_res honors single-model selection parameters", {
  skip_if_not_installed("glmnet")

  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "M")
  candidate_genes <- paste0("G", seq_len(6))
  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  processed_train <- iklSurvML:::preprocess_train_data(train_data, common_feature)
  est_dd <- as.data.frame(processed_train)[, common_feature[-1]]

  lasso_fit <- iklSurvML:::train_lasso(est_dd, candidate_genes, seed = 5201314)
  ridge_fit <- iklSurvML:::train_ridge(est_dd, candidate_genes, seed = 5201314)
  result <- list(
    ml.res = list(Lasso = lasso_fit, Ridge = ridge_fit),
    Sig.genes = candidate_genes
  )

  rs <- iklSurvML::cal_RS_ml_res(
    res.by.ML.Dev.Prog.Sig = result,
    train_data = train_data,
    inputmatrix.list = list(Train = train_data),
    mode = "single",
    single_ml = "Ridge"
  )

  expect_equal(names(rs), "Ridge")
})

test_that("survival model results carry explicit prediction metadata", {
  skip_if_not_installed("glmnet")

  train_data <- make_survival_smoke_data(n = 35, p = 6, prefix = "MI")
  candidate_genes <- paste0("G", seq_len(6))

  result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
    train_data = train_data,
    list_train_vali_Data = list(Train = train_data),
    candidate_genes = candidate_genes,
    unicox.filter.for.candi = FALSE,
    mode = "single",
    single_ml = "Lasso",
    use_parallel = FALSE
  )

  expect_true("Model.info" %in% names(result))
  expect_equal(result$Model.info$Lasso$pipeline$learner, "Lasso")
  expect_equal(result$Model.info$Lasso$features, candidate_genes)
  expect_equal(result$Model.info$Lasso$risk_direction, "higher_is_worse")
})

test_that("legacy double-mode RSF second-stage uses caller nodesize", {
  combo_bodies <- paste(
    deparse(iklSurvML:::run_stepcox_combination),
    deparse(iklSurvML:::run_coxboost_combination),
    deparse(iklSurvML:::run_lasso_combination),
    collapse = "\n"
  )

  expect_false(grepl("train_rsf\\(est_dd2, 5, seed\\)", combo_bodies))
  expect_true(grepl("train_rsf\\(est_dd2, rf_nodesize, seed\\)", combo_bodies))
})

test_that("all-mode implementations expose the fixed 117 grid", {
  expect_true(exists("run_all_algorithms_128", envir = asNamespace("iklSurvML"), inherits = FALSE))
  expect_true(exists("run_all_algorithms_128_parallel", envir = asNamespace("iklSurvML"), inherits = FALSE))

  legacy_body <- paste(deparse(iklSurvML:::run_all_algorithms), collapse = "\n")
  optimized_body <- paste(deparse(iklSurvML:::run_all_algorithms_128), collapse = "\n")
  parallel_body <- paste(deparse(iklSurvML:::run_all_algorithms_128_parallel), collapse = "\n")

  expect_equal(iklSurvML:::normalize_all_mode_model_grid(NULL), "117")
  expect_equal(iklSurvML:::normalize_all_mode_model_grid("117"), "117")
  expect_equal(iklSurvML:::all_mode_model_grid_size(), 117L)
  expect_equal(iklSurvML:::all_mode_model_grid_size("117"), 117L)
  expect_equal(iklSurvML:::all_mode_stepcox_selector_dirs(), c("both", "backward", "forward"))
  expect_false("RSF" %in% iklSurvML:::all_mode_coxboost_second_stage_algorithms())
  expect_false(any(c("Enet", "Ridge") %in% iklSurvML:::all_mode_lasso_second_stage_algorithms()))
  expect_error(iklSurvML:::normalize_all_mode_model_grid("101"), "only supports '117'")
  expect_error(iklSurvML:::normalize_all_mode_model_grid("128"), "only supports '117'")

  expect_true(grepl("all_mode_model_grid_size", legacy_body, fixed = TRUE))
  expect_true(grepl("all_mode_stepcox_selector_dirs", legacy_body, fixed = TRUE))

  for (body_text in list(optimized_body, parallel_body)) {
    expect_true(grepl("all_mode_model_grid_size", body_text, fixed = TRUE))
    expect_true(grepl("all_mode_stepcox_selector_dirs", body_text, fixed = TRUE))
    expect_false(grepl("CoxBoost + RSF", body_text, fixed = TRUE))
    expect_false(grepl("Lasso + Enet", body_text, fixed = TRUE))
    expect_false(grepl("Lasso + Ridge", body_text, fixed = TRUE))
    expect_false(grepl("RSF + RSF", body_text, fixed = TRUE))
    expect_false(grepl("CoxBoost + CoxBoost", body_text, fixed = TRUE))
    expect_false(grepl("Lasso + Lasso", body_text, fixed = TRUE))
    expect_false(grepl("StepCox\\[[^]]+\\] \\+ StepCox", body_text))
  }
})

test_that("survival algorithm validation reads from the central registry", {
  registry <- iklSurvML:::survival_ml_registry()

  expect_equal(iklSurvML:::survival_ml_names(), registry$name)
  expect_equal(
    iklSurvML:::survival_ml_first_stage_names(),
    c("RSF", "StepCox", "CoxBoost", "Lasso")
  )
  expect_true(all(c("name", "display", "first_stage_selector") %in% colnames(registry)))
})

test_that("previous prognostic signature helpers reject missing genes instead of zero-imputing them", {
  sig <- data.frame(
    model = "MockSig",
    PMID = "0",
    Cancer = "Mock",
    Author = "Mock",
    Coef = 1,
    symbol = "MissingGene",
    stringsAsFactors = FALSE
  )
  cohort <- data.frame(
    ID = paste0("S", seq_len(5)),
    OS.time = seq(100, 500, by = 100),
    OS = c(1, 0, 1, 0, 1),
    PresentGene = seq_len(5),
    check.names = FALSE
  )

  expect_error(
    iklSurvML::cal_RS_pre.prog.sig(
      use_your_own_collected_sig = TRUE,
      collected_sig_table = sig,
      type.sig = NULL,
      list_input_data = list(Cohort = cohort)
    ),
    "not imputed as zero"
  )
  expect_error(
    iklSurvML::cal_cindex_pre.prog.sig(
      use_your_own_collected_sig = TRUE,
      collected_sig_table = sig,
      type.sig = NULL,
      list_input_data = list(Cohort = cohort)
    ),
    "not imputed as zero"
  )
  expect_error(
    iklSurvML::cal_auc_pre.prog.sig(
      use_your_own_collected_sig = TRUE,
      collected_sig_table = sig,
      type.sig = NULL,
      list_input_data = list(Cohort = cohort),
      AUC_time = 1
    ),
    "not imputed as zero"
  )
})

test_that("cal_AUC_ml_res fails clearly when AUC_time exceeds follow-up", {
  train_data <- make_survival_smoke_data(n = 10, p = 2, prefix = "A")
  train_data$OS.time <- seq(10, 100, by = 10)
  result <- list(Sig.genes = c("G1", "G2"), ml.res = list())

  expect_error(
    iklSurvML::cal_AUC_ml_res(
      res.by.ML.Dev.Prog.Sig = result,
      train_data = train_data,
      inputmatrix.list = list(Train = train_data),
      mode = "all",
      AUC_time = 1
    ),
    "outside the available follow-up window"
  )
})

test_that("signature_score PCA branch uses the caller signature and removes constant genes", {
  suppressWarnings(skip_if_not_installed("IOBR"))

  cohort <- data.frame(
    ID = paste0("S", seq_len(6)),
    OS.time = seq(100, 600, by = 100),
    OS = c(1, 0, 1, 0, 1, 0),
    ConstantGene = rep(1, 6),
    G1 = seq_len(6),
    G2 = c(2, 3, 5, 7, 11, 13),
    check.names = FALSE
  )
  result <- list(
    riskscore = list(
      MockModel = list(
        Cohort = data.frame(
          ID = cohort$ID,
          OS.time = cohort$OS.time,
          OS = cohort$OS,
          RS = seq_len(nrow(cohort))
        )
      )
    )
  )

  out <- NULL
  expect_no_error({
    out <- suppressWarnings(
      iklSurvML::signature_score(
        res.by.ML.Dev.Prog.Sig = result,
        inputmatrix.list = list(Cohort = cohort),
        signature = list(MockSignature = c("G1", "G2")),
        estima_method = "pca",
        diff.test = FALSE,
        correlation.test = FALSE
      )
    )
  })
  expect_equal(names(out), c("signature_list", "differential_sig_list"))
})

test_that("ML.Corefeature.Prog.Screen fails clearly when the mandatory prefilter selects no genes", {
  input <- data.frame(
    ID = paste0("S", seq_len(8)),
    OS.time = seq_len(8) + 10,
    OS = rep(c(1, 0), 4),
    G1 = rep(1, 8),
    G2 = rep(2, 8),
    check.names = FALSE
  )

  expect_error(
    iklSurvML::ML.Corefeature.Prog.Screen(
      InputMatrix = input,
      candidate_genes = c("G1", "G2"),
      mode = "single",
      single_ml = "StepCox",
      seed = 1,
      nodesize = 5
    ),
    "No genes passed the univariate Cox prefilter"
  )
})

test_that("TME_deconvolution_all handles vector microarray_names without scalar-if logic", {
  body_text <- paste(deparse(iklSurvML::TME_deconvolution_all), collapse = "\n")

  expect_false(grepl('if \\(selected_columns == "none"\\)', body_text))
  expect_true(grepl("microarray_is_none", body_text, fixed = TRUE))
  expect_true(grepl("invalid_microarray_names", body_text, fixed = TRUE))
})

test_that("TME_deconvolution_all reports cohort failures instead of returning silent NULLs", {
  cohort <- data.frame(
    ID = paste0("S", seq_len(4)),
    OS.time = seq(100, 400, by = 100),
    OS = c(1, 0, 1, 0),
    GeneA = c(1, 2, 3, 4),
    GeneB = c(2, 3, 4, 5),
    check.names = FALSE
  )

  expect_error(
    iklSurvML::TME_deconvolution_all(
      inputmatrix.list = list(Cohort = cohort),
      deconvolution_method = "not_a_method"
    ),
    "TME deconvolution failed"
  )
})

test_that("backward compatibility aliases are exported", {
  exports <- getNamespaceExports("iklSurvML")
  expect_true(all(c(
    "cal_auc_ml_res",
    "cal_rs_ml_res",
    "cal_rs_pre_prog_sig",
    "cal_auc_pre_prog_sig",
    "cal_cindex_pre_prog_sig",
    "cal_auc_category_model",
    "cal_roc_category_model"
  ) %in% exports))
})

test_that("Xgboost screening has an explicit optional dependency gate", {
  body_text <- paste(deparse(iklSurvML::ML.Corefeature.Prog.Screen), collapse = "\n")

  expect_true(grepl('requireNamespace\\("xgboost"', body_text))
  expect_true(grepl("install xgboost", body_text, fixed = TRUE))
})

test_that("core-feature survival screening reuses central survival model wrappers", {
  body_text <- paste(deparse(iklSurvML::ML.Corefeature.Prog.Screen), collapse = "\n")

  expect_true(grepl("train_rsf", body_text, fixed = TRUE))
  expect_true(grepl("get_rsf_selected_vars", body_text, fixed = TRUE))
  expect_true(grepl("train_coxboost", body_text, fixed = TRUE))
  expect_true(grepl("get_coxboost_selected_vars", body_text, fixed = TRUE))
  expect_true(grepl("train_stepcox", body_text, fixed = TRUE))
  expect_true(grepl("get_stepcox_selected_vars", body_text, fixed = TRUE))
  expect_false(grepl("optimCoxBoostPenalty", body_text, fixed = TRUE))
  expect_false(grepl("stats::step\\(coxph\\(Surv\\(OS.time, OS\\) ~ \\\\., est_dd\\)", body_text))
})

test_that("core-feature screening avoids unsafe zero imputation and removes Boruta", {
  body_text <- paste(deparse(iklSurvML::ML.Corefeature.Prog.Screen), collapse = "\n")

  expect_true(grepl("fit_survival_preprocess_recipe", body_text, fixed = TRUE))
  expect_false(grepl("inputSet[is.na(inputSet)] <- 0", body_text, fixed = TRUE))
  expect_false(grepl("y = survival::Surv", body_text, fixed = TRUE))
  expect_false(grepl("screen_boruta_survival_vars", body_text, fixed = TRUE))
  expect_false(grepl("single_ml == \"Boruta\"", body_text, fixed = TRUE))
  expect_false(grepl("append_screen_result(selected.feature, \"Boruta\"", body_text, fixed = TRUE))
})

test_that("core-feature screening rejects removed Boruta selector before fitting", {
  train_data <- make_survival_smoke_data(n = 12, p = 3, prefix = "B")

  expect_error(
    iklSurvML::ML.Corefeature.Prog.Screen(
      InputMatrix = train_data,
      candidate_genes = paste0("G", seq_len(3)),
      mode = "single",
      seed = 5201314,
      single_ml = "Boruta"
    ),
    "Boruta"
  )
})

test_that("immunotherapy signature scoring validates named genes and avoids positional slices", {
  cohort <- data.frame(
    ID = paste0("S", seq_len(3)),
    Var = c("Y", "N", "Y"),
    G2 = c(2, 4, 6),
    G1 = c(1, 3, 5),
    check.names = FALSE
  )

  expect_equal(
    iklSurvML:::calculate_signature_score_by_genes(
      cohort,
      signature_genes = c("G1", "G2"),
      signature_name = "Mock.Sig",
      dataset_name = "MockCohort",
      score_type = "mean"
    ),
    rowMeans(cohort[, c("G1", "G2")])
  )
  expect_error(
    iklSurvML:::calculate_signature_score_by_genes(
      cohort,
      signature_genes = c("G1", "MissingGene"),
      signature_name = "Mock.Sig",
      dataset_name = "MockCohort",
      score_type = "mean"
    ),
    "missing required genes"
  )

  body_text <- paste(deparse(iklSurvML::cal_auc_previous_sig), collapse = "\n")
  expect_false(grepl("rowMeans\\(new\\[, 2:7\\]\\)", body_text))
  expect_false(grepl("new\\$PDCD1", body_text, fixed = TRUE))
})

test_that("all-mode count shortfalls warn instead of stopping", {
  expect_warning(
    expect_invisible(iklSurvML:::warn_if_all_mode_incomplete(116L, expected = 117L, context = "test all-mode")),
    "produced 116 models; expected 117"
  )
  expect_silent(iklSurvML:::warn_if_all_mode_incomplete(117L, expected = 117L, context = "test all-mode"))
})

test_that("all-mode partial results require an explicit opt-in", {
  partial_result <- list(
    ml.res = list(Mock = structure(list(), class = "mock")),
    Model.errors = "Mock + GBM: backend failed"
  )

  expect_error(
    iklSurvML:::assert_complete_all_mode_result(
      partial_result,
      expected = 117L,
      allow_partial = FALSE,
      context = "test all-mode"
    ),
    "allow_partial = TRUE"
  )
  expect_silent(
    iklSurvML:::assert_complete_all_mode_result(
      partial_result,
      expected = 117L,
      allow_partial = TRUE,
      context = "test all-mode"
    )
  )
})
