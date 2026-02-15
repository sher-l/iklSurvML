# Machine Learning Algorithms for Prognosis
#
# This file contains core implementations of 10 ML algorithms:
# RSF, Enet, StepCox, CoxBoost, plsRcox, superpc, GBM, survivalsvm, Ridge, Lasso

# ---- RSF (Random Survival Forest) ----

#' Train RSF model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param rf_nodesize Node size parameter (default 5)
#' @param seed Random seed
#' @return Trained RSF model
#' @keywords internal
train_rsf <- function(est_dd, rf_nodesize = 5, seed = 5201314) {
  set.seed(seed)
  # Use survival::Surv explicitly and construct formula in calling environment
  Surv <- survival::Surv
  fit <- randomForestSRC::rfsrc(
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
  return(fit)
}

#' Get RSF selected variables
#'
#' @param fit RSF model
#' @return Character vector of selected variable names
#' @keywords internal
get_rsf_selected_vars <- function(fit) {
  rid <- randomForestSRC::var.select(object = fit, conservative = "high")
  return(rid$topvars)
}

#' Predict with RSF model
#'
#' @param fit RSF model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_rsf <- function(fit, newdata) {
  return(predict(fit, newdata = newdata)$predicted)
}

# ---- Enet (Elastic Net) ----

#' Train Enet model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param rid Feature names to use (optional)
#' @param alpha Alpha parameter for elastic net (0-1)
#' @param seed Random seed
#' @return Trained cv.glmnet model
#' @keywords internal
train_enet <- function(est_dd, rid = NULL, alpha = 0.5, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  fit <- glmnet::cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  return(fit)
}

#' Predict with Enet model
#'
#' @param fit cv.glmnet model
#' @param newdata New data for prediction
#' @param rid Feature names to use
#' @return Predicted risk scores
#' @keywords internal
predict_enet <- function(fit, newdata, rid) {
  return(as.numeric(predict(
    fit,
    type = "link",
    newx = as.matrix(newdata[, rid]),
    s = fit$lambda.min
  )))
}

# ---- Lasso ----

#' Train Lasso model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param rid Feature names to use (optional)
#' @param seed Random seed
#' @return Trained cv.glmnet model
#' @keywords internal
train_lasso <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x1, x2,
    nfold = 10,
    family = "cox",
    alpha = 1
  )
  return(fit)
}

#' Get Lasso selected variables
#'
#' @param fit cv.glmnet model
#' @return Character vector of selected variable names
#' @keywords internal
get_lasso_selected_vars <- function(fit) {
  my_coefs <- coef(fit, s = "lambda.min")
  rid <- my_coefs@Dimnames[[1]][Matrix::which(my_coefs != 0)]
  return(rid)
}

#' Predict with Lasso model
#'
#' @param fit cv.glmnet model
#' @param newdata New data for prediction
#' @param rid Feature names to use
#' @return Predicted risk scores
#' @keywords internal
predict_lasso <- function(fit, newdata, rid) {
  return(as.numeric(predict(
    fit,
    type = "response",
    newx = as.matrix(newdata[, rid]),
    s = fit$lambda.min
  )))
}

# ---- Ridge ----

#' Train Ridge model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param rid Feature names to use (optional)
#' @param seed Random seed
#' @return Trained cv.glmnet model
#' @keywords internal
train_ridge <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  # Match original code: call glmnet first (consumes random numbers), then cv.glmnet
  fit <- glmnet::glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
  cv_fit <- glmnet::cv.glmnet(x1, x2, nfold = 10, family = "cox")
  return(list(fit = fit, cv.fit = cv_fit))
}

#' Predict with Ridge model
#'
#' @param fit cv.glmnet model
#' @param newdata New data for prediction
#' @param rid Feature names to use
#' @return Predicted risk scores
#' @keywords internal
predict_ridge <- function(fit, newdata, rid) {
  # fit is a list with $fit and $cv.fit (matching original code structure)
  # Original code: predict(fit, s = cv.fit$lambda.min) - uses glmnet object, not cv.glmnet
  if (is.list(fit) && "fit" %in% names(fit)) {
    glmnet_fit <- fit$fit
    lambda <- fit$cv.fit$lambda.min
  } else {
    glmnet_fit <- fit
    lambda <- fit$lambda.min
  }
  return(as.numeric(predict(
    glmnet_fit,
    type = "response",
    newx = as.matrix(newdata[, rid]),
    s = lambda
  )))
}

# ---- StepCox ----

#' Train StepCox model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param direction Direction for stepwise selection: "both", "backward", "forward"
#' @return Trained coxph model
#' @keywords internal
train_stepcox <- function(est_dd, direction = "both") {
  fit <- step(
    survival::coxph(survival::Surv(OS.time, OS) ~ ., est_dd),
    direction = direction
  )
  return(fit)
}

#' Get StepCox selected variables
#'
#' @param fit coxph model from stepcox
#' @return Character vector of selected variable names
#' @keywords internal
get_stepcox_selected_vars <- function(fit) {
  return(names(coef(fit)))
}

#' Predict with StepCox model
#'
#' @param fit coxph model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_stepcox <- function(fit, newdata) {
  return(predict(fit, type = "risk", newdata = newdata))
}

# ---- CoxBoost ----

#' Train CoxBoost model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @return Trained CoxBoost model
#' @keywords internal
train_coxboost <- function(est_dd, seed = 5201314) {
  set.seed(seed)
  pen <- CoxBoost::optimCoxBoostPenalty(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    trace = TRUE,
    start.penalty = 500,
    parallel = TRUE
  )

  cv_res <- CoxBoost::cv.CoxBoost(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    maxstepno = 500,
    K = 10,
    type = "verweij",
    penalty = pen$penalty
  )

  fit <- CoxBoost::CoxBoost(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    stepno = cv_res$optimal.step,
    penalty = pen$penalty
  )
  return(fit)
}

#' Predict with CoxBoost model
#'
#' @param fit CoxBoost model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_coxboost <- function(fit, newdata) {
  return(as.numeric(predict(
    fit,
    newdata = newdata[, -c(1, 2)],
    newtime = newdata[, 1],
    newstatus = newdata[, 2],
    type = "lp"
  )))
}

# ---- plsRcox ----

#' Train plsRcox model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param rid Feature names to use (optional)
#' @param seed Random seed
#' @return Trained plsRcox model
#' @keywords internal
train_plsrcox <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  set.seed(seed)
  cv_plsrcox_res <- plsRcox::cv.plsRcox(
    list(
      x = est_dd[, rid],
      time = est_dd$OS.time,
      status = est_dd$OS
    ),
    nt = 10,
    verbose = FALSE
  )

  fit <- plsRcox::plsRcox(
    est_dd[, rid],
    time = est_dd$OS.time,
    event = est_dd$OS,
    nt = as.numeric(cv_plsrcox_res[5])
  )
  return(fit)
}

#' Predict with plsRcox model
#'
#' @param fit plsRcox model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_plsrcox <- function(fit, newdata) {
  return(as.numeric(predict(
    fit,
    type = "lp",
    newdata = newdata[, -c(1, 2)]
  )))
}

# ---- SuperPC ----

#' Train SuperPC model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @return List with fit and cv.fit
#' @keywords internal
train_superpc <- function(est_dd, seed = 5201314) {
  data <- list(
    x = t(est_dd[, -c(1, 2)]),
    y = est_dd$OS.time,
    censoring.status = est_dd$OS,
    featurenames = colnames(est_dd)[-c(1, 2)]
  )
  set.seed(seed)
  fit <- superpc::superpc.train(
    data = data,
    type = "survival",
    s0.perc = 0.5
  )

  # Retry on error
  repeat {
    tryCatch({
      cv_fit <- superpc::superpc.cv(
        fit, data,
        n.threshold = 20,
        n.fold = 10,
        n.components = 3,
        min.features = 2,
        max.features = nrow(data$x),
        compute.fullcv = TRUE,
        compute.preval = TRUE
      )
      break
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      cat("Retrying...\n")
      Sys.sleep(1)
    })
  }

  return(list(fit = fit, cv_fit = cv_fit))
}

#' Predict with SuperPC model
#'
#' @param fit SuperPC model
#' @param cv_fit Cross-validation result
#' @param train_data Training data used for model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_superpc <- function(fit, cv_fit, train_data, newdata) {
  data <- list(
    x = t(train_data[, -c(1, 2)]),
    y = train_data$OS.time,
    censoring.status = train_data$OS,
    featurenames = colnames(train_data)[-c(1, 2)]
  )

  test <- list(
    x = t(newdata[, -c(1, 2)]),
    y = newdata$OS.time,
    censoring.status = newdata$OS,
    featurenames = colnames(newdata)[-c(1, 2)]
  )

  ff <- superpc::superpc.predict(
    fit, data, test,
    threshold = cv_fit$thresholds[which.max(cv_fit[["scor"]][1, ])],
    n.components = 1
  )
  return(as.numeric(ff$v.pred))
}

# ---- GBM ----

#' Train GBM model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @param cores_for_parallel Number of cores for parallel processing
#' @return List with fit and best iteration
#' @keywords internal
train_gbm <- function(est_dd, seed = 5201314, cores_for_parallel = 6) {
  set.seed(seed)
  fit <- gbm::gbm(
    formula = survival::Surv(OS.time, OS) ~ .,
    data = est_dd,
    distribution = "coxph",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = cores_for_parallel
  )

  # Find optimal number of trees
  best <- which.min(fit$cv.error)

  set.seed(seed)
  fit <- gbm::gbm(
    formula = survival::Surv(OS.time, OS) ~ .,
    data = est_dd,
    distribution = "coxph",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = cores_for_parallel
  )

  return(list(fit = fit, best = best))
}

#' Predict with GBM model
#'
#' @param fit GBM model
#' @param best Best iteration number
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_gbm <- function(fit, best, newdata) {
  return(as.numeric(predict(fit, newdata, n.trees = best, type = "link")))
}

# ---- survivalsvm ----

#' Train survivalsvm model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @return Trained survivalsvm model
#' @keywords internal
train_survivalsvm <- function(est_dd) {
  fit <- survivalsvm::survivalsvm(
    survival::Surv(OS.time, OS) ~ .,
    data = est_dd,
    gamma.mu = 1
  )
  return(fit)
}

#' Predict with survivalsvm model
#'
#' @param fit survivalsvm model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
predict_survivalsvm <- function(fit, newdata) {
  return(as.numeric(predict(fit, newdata)$predicted))
}

# ---- Utility functions ----

#' Calculate C-index for predictions
#'
#' @param rs Risk scores
#' @param data Data with OS.time and OS
#' @return C-index value
#' @keywords internal
calculate_cindex <- function(rs, data) {
  tmp <- data.frame(RS = rs, OS.time = data$OS.time, OS = data$OS)
  return(as.numeric(summary(
    survival::coxph(survival::Surv(OS.time, OS) ~ RS, tmp)
  )$concordance[1]))
}

#' Calculate risk scores for validation data
#'
#' @param val_dd_list List of validation data
#' @param rs_func Function to calculate risk scores
#' @return List of risk score data frames
#' @keywords internal
calculate_risk_scores <- function(val_dd_list, rs_func) {
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = rs_func(x))
  })
  return(rs)
}
