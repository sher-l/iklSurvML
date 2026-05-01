test_that("category ML default methods exclude unavailable optional engines", {
  methods <- iklSurvML:::resolve_category_methods(NULL)

  if (!requireNamespace("cancerclass", quietly = TRUE)) {
    expect_false("cancerclass" %in% methods)
  }
  if (!requireNamespace("klaR", quietly = TRUE)) {
    expect_false("nb" %in% methods)
  }
  if (!requireNamespace("fastAdaboost", quietly = TRUE)) {
    expect_false("adaboost" %in% methods)
  }
})

test_that("category method dependency map includes transitive cancerclass runtime packages", {
  deps <- iklSurvML:::category_method_packages()

  expect_true(all(c("cancerclass", "Biobase", "pROC") %in% deps$cancerclass))
})

test_that("requesting an unavailable optional category engine fails clearly", {
  skip_if(requireNamespace("cancerclass", quietly = TRUE), "cancerclass is installed")

  expect_error(
    iklSurvML:::resolve_category_methods("cancerclass"),
    "cancerclass"
  )
})

make_category_smoke_data <- function(n = 45, p = 6, prefix = "C") {
  set.seed(n + p)
  genes <- paste0("G", seq_len(p))
  expr <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(expr) <- genes
  linpred <- expr[, 1] - 0.8 * expr[, 2]
  data.frame(
    ID = paste0(prefix, seq_len(n)),
    Var = ifelse(linpred + stats::rnorm(n, sd = 0.5) > 0, "Y", "N"),
    expr,
    check.names = FALSE
  )
}

test_that("category ML imputes training feature NAs from a train-only recipe", {
  skip_if_not_installed("randomForest")

  train_data <- make_category_smoke_data(n = 45, p = 6, prefix = "T")
  validation_data <- make_category_smoke_data(n = 25, p = 6, prefix = "V")
  train_data$G1[3] <- NA_real_
  validation_data$G1[4] <- NA_real_

  warnings_seen <- character()
  result <- withCallingHandlers(
    iklSurvML::ML.Dev.Pred.Category.Sig(
      train_data = train_data,
      list_train_vali_Data = list(validation = validation_data),
      candidate_genes = paste0("G", seq_len(6)),
      methods = "rf",
      seed = 123,
      cores_for_parallel = 1
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(names(result$model), "rf")
  expect_s3_class(result$Preprocess.recipe, "ikl_category_preprocess_recipe")
  expect_equal(
    result$Preprocess.recipe$impute_values[["G1"]],
    mean(train_data$G1, na.rm = TRUE)
  )
  expect_false(any(grepl("missing values in object", warnings_seen, fixed = TRUE)))
})

test_that("category preprocessing rejects duplicated sample IDs", {
  train_data <- make_category_smoke_data(n = 12, p = 4, prefix = "DUP")
  train_data$ID[2] <- train_data$ID[1]

  expect_error(
    iklSurvML:::fit_category_preprocess_recipe(
      train_data,
      common_feature = c("ID", "Var", paste0("G", seq_len(4))),
      positive_class = "Y"
    ),
    "duplicated sample IDs"
  )
})

test_that("category ML uses Y as the default positive class", {
  skip_if_not_installed("randomForest")

  train_data <- make_category_smoke_data(n = 45, p = 6, prefix = "T")
  validation_data <- make_category_smoke_data(n = 25, p = 6, prefix = "V")

  result <- iklSurvML::ML.Dev.Pred.Category.Sig(
    train_data = train_data,
    list_train_vali_Data = list(validation = validation_data),
    candidate_genes = paste0("G", seq_len(6)),
    methods = "rf",
    seed = 123,
    cores_for_parallel = 1
  )

  expect_equal(result$positive_class, "Y")
  expect_identical(as.character(result$model$rf$levels), c("Y", "N"))
})

test_that("category RF tuning grid stays within the available feature count", {
  grid <- iklSurvML:::category_tune_grid(n_features = 6)

  expect_true(all(grid$rf$mtry >= 1))
  expect_true(all(grid$rf$mtry <= 6))
  expect_equal(length(unique(grid$rf$mtry)), length(grid$rf$mtry))
})

test_that("requesting an unavailable caret category engine fails clearly", {
  skip_if(requireNamespace("klaR", quietly = TRUE), "klaR is installed")

  expect_error(
    iklSurvML:::resolve_category_methods("nb"),
    "klaR"
  )
})

test_that("category ML input validation fails before starting workers", {
  invalid_train <- data.frame(
    ID = c("S1", "S2"),
    Outcome = c("Y", "N"),
    G1 = c(1, 2),
    G2 = c(2, 1),
    check.names = FALSE
  )

  expect_error(
    iklSurvML::ML.Dev.Pred.Category.Sig(
      train_data = invalid_train,
      list_train_vali_Data = list(Train = invalid_train),
      candidate_genes = c("G1", "G2"),
      methods = "nb",
      cores_for_parallel = 1
    ),
    "feature_alignment='strict'|first 2 columns"
  )
})

test_that("category ML strict feature alignment rejects missing validation features", {
  train_data <- make_category_smoke_data(n = 20, p = 4, prefix = "T")
  validation_data <- make_category_smoke_data(n = 15, p = 4, prefix = "V")
  validation_data$G4 <- NULL

  expect_error(
    iklSurvML::ML.Dev.Pred.Category.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Val = validation_data),
      candidate_genes = paste0("G", seq_len(4)),
      methods = "rf",
      cores_for_parallel = 1
    ),
    "feature_alignment='strict'"
  )

  skip_if_not_installed("randomForest")
  expect_warning(
    result <- iklSurvML::ML.Dev.Pred.Category.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Val = validation_data),
      candidate_genes = paste0("G", seq_len(4)),
      methods = "rf",
      cores_for_parallel = 1,
      feature_alignment = "intersection"
    ),
    "dropped candidate features"
  )
  expect_equal(result$sig.gene, paste0("G", seq_len(3)))
})

test_that("category ML validates worker count before creating a PSOCK cluster", {
  train_data <- make_category_smoke_data(n = 20, p = 4, prefix = "T")

  expect_error(
    iklSurvML::ML.Dev.Pred.Category.Sig(
      train_data = train_data,
      list_train_vali_Data = list(Train = train_data),
      candidate_genes = paste0("G", seq_len(4)),
      methods = "rf",
      cores_for_parallel = 0
    ),
    "cores_for_parallel"
  )
})

test_that("exported category AUC/ROC helpers accept caret train models", {
  skip_if_not_installed("randomForest")

  train_data <- make_category_smoke_data(n = 45, p = 6, prefix = "TA")
  validation_data <- make_category_smoke_data(n = 25, p = 6, prefix = "VA")

  result <- iklSurvML::ML.Dev.Pred.Category.Sig(
    train_data = train_data,
    list_train_vali_Data = list(validation = validation_data),
    candidate_genes = paste0("G", seq_len(6)),
    methods = "rf",
    seed = 123,
    cores_for_parallel = 1
  )

  expect_no_error(auc <- iklSurvML::cal.auc.category.model(result, validation_data))
  expect_no_error(roc <- iklSurvML::cal.roc.category.model(result, validation_data))
  expect_equal(rownames(auc), "rf")
  expect_equal(names(roc), "rf")
})

test_that("category tuning defaults avoid the exhaustive SVM grid", {
  standard <- iklSurvML:::category_tune_grid(n_features = 6)
  exhaustive <- iklSurvML:::category_tune_grid(n_features = 6, tune_profile = "exhaustive")

  expect_lt(nrow(standard$svmRadialWeights), nrow(exhaustive$svmRadialWeights))
  expect_lte(nrow(standard$svmRadialWeights), 27)
  expect_equal(nrow(exhaustive$svmRadialWeights), 175)
})

test_that("exported category AUC/ROC helpers fail clearly for missing features", {
  result <- list(
    model = list(rf = structure(list(), class = "mock")),
    sig.gene = c("G1", "MissingGene"),
    positive_class = "Y"
  )
  cohort <- make_category_smoke_data(n = 10, p = 2, prefix = "MISS")

  expect_error(
    iklSurvML::cal.auc.category.model(result, cohort),
    "missing required columns"
  )
  expect_error(
    iklSurvML::cal.roc.category.model(result, cohort),
    "missing required columns"
  )
})
