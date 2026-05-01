make_b2_smoke_fixture <- function(seed = 101, n_train = 50, n_val = 30, p = 20) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(p))

  make_matrix <- function(n, prefix) {
    signal <- stats::rnorm(n)
    expr <- vapply(seq_len(p), function(j) {
      signal * cos(j / 4) + stats::rnorm(n, sd = 0.45 + j / 100)
    }, numeric(n))
    colnames(expr) <- genes
    linpred <- 0.8 * expr[, 1] - 0.6 * expr[, 2] + 0.5 * expr[, 3]
    data.frame(
      ID = paste0(prefix, seq_len(n)),
      OS.time = pmax(30, round(exp(8 - linpred + stats::rnorm(n, sd = 0.35)), 3)),
      OS = as.integer(stats::runif(n) < stats::plogis(linpred)),
      expr,
      check.names = FALSE
    )
  }

  list(
    train_data = make_matrix(n_train, "S"),
    list_train_vali_Data = list(validation = make_matrix(n_val, "V")),
    candidate_genes = genes,
    seed = 5201314,
    nodesize = 5,
    meta = list(samples_train = n_train, samples_val = n_val, genes = p, generator_seed = seed)
  )
}

shared_smoke_path <- Sys.getenv("IKLSURVML_B2_FIXTURE", unset = NA_character_)

load_shared_smoke_inputs <- function() {
  dataset <- if (!is.na(shared_smoke_path) && nzchar(shared_smoke_path) && file.exists(shared_smoke_path)) {
    readRDS(shared_smoke_path)
  } else {
    make_b2_smoke_fixture()
  }

  train_data <- dataset$train_data
  list_train_vali_Data <- dataset$list_train_vali_Data
  candidate_genes <- dataset$candidate_genes

  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  list(
    dataset = dataset,
    train_data = iklSurvML:::preprocess_train_data(train_data, common_feature),
    list_train_vali_Data = iklSurvML:::preprocess_data_list(list_train_vali_Data, common_feature),
    candidate_genes = candidate_genes,
    common_feature = common_feature
  )
}

test_that("train_survivalsvm falls back for the shared Lasso subspace", {
  inputs <- load_shared_smoke_inputs()
  est_dd <- as.data.frame(inputs$train_data)[, inputs$common_feature[-1]]
  pre_var <- inputs$common_feature[-c(1:3)]

  lasso_fit <- iklSurvML:::train_lasso(est_dd, pre_var, inputs$dataset$seed)
  lasso_vars <- iklSurvML:::get_lasso_selected_vars(lasso_fit)
  est_lasso <- inputs$train_data[, c("OS.time", "OS", lasso_vars)]

  fit <- iklSurvML:::train_survivalsvm(est_lasso, inputs$dataset$seed)

  expect_s3_class(fit, "survivalsvm")
})

test_that("parallel all-mode keeps B2 materialized models and reaches the fixed 117-grid", {
  skip_on_os("windows")
  skip_if_not(
    identical(Sys.getenv("RUN_LONG_B2_REGRESSION"), "true"),
    "long-running integration regression covered by scripts/run_b2_retest_12c.R"
  )
  skip_if(
    is.na(shared_smoke_path) || !nzchar(shared_smoke_path) || !file.exists(shared_smoke_path),
    "set IKLSURVML_B2_FIXTURE to run the full 117-model B2 regression"
  )
  inputs <- load_shared_smoke_inputs()

  result <- NULL
  invisible(capture.output({
    result <- iklSurvML::ML.Dev.Prog.Sig.Fast(
      train_data = inputs$dataset$train_data,
      list_train_vali_Data = inputs$dataset$list_train_vali_Data,
      candidate_genes = inputs$dataset$candidate_genes,
      unicox.filter.for.candi = FALSE,
      mode = "all",
      seed = inputs$dataset$seed,
      nodesize = inputs$dataset$nodesize,
      use_parallel = TRUE,
      cores_for_parallel = 12
    )
  }))

  model_names <- names(result$ml.res)
  expect_equal(length(model_names), 117)
  expect_true(all(c(
    "StepCox[both] + RSF",
    "StepCox[backward] + RSF",
    "StepCox[forward] + RSF",
    "RSF + CoxBoost",
    "Lasso + survival-SVM"
  ) %in% model_names))
  expect_false("CoxBoost + RSF" %in% model_names)
  expect_false("Lasso + Enet[\u03b1=0.1]" %in% model_names)
  expect_false("Lasso + Ridge" %in% model_names)
})
