shared_smoke_path <- "/tmp/iklSurvML-audit/shared-smoke-50x20-seed101.rds"

load_shared_smoke_inputs <- function() {
  skip_if_not(file.exists(shared_smoke_path), paste("missing fixture:", shared_smoke_path))
  dataset <- readRDS(shared_smoke_path)

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

test_that("parallel all-mode keeps B2 materialized models and reaches 117", {
  skip_on_os("windows")
  skip_if_not(
    identical(Sys.getenv("RUN_LONG_B2_REGRESSION"), "true"),
    "long-running integration regression covered by scripts/run_b2_retest_12c.R"
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
    "Lasso + survival-SVM"
  ) %in% model_names))
})
