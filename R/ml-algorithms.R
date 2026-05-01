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
#' @noRd
train_rsf <- function(est_dd, rf_nodesize = 5, seed = 5201314) {
  set.seed(seed)
  features <- colnames(est_dd)[-c(1, 2)]
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
  attr(fit, "iklsurvml_features") <- features
  return(fit)
}

#' Get RSF selected variables
#'
#' @param fit RSF model
#' @param seed Random seed for variable selection
#' @return Character vector of selected variable names
#' @keywords internal
#' @noRd
get_rsf_selected_vars <- function(fit, seed = 5201314) {
  set.seed(seed)
  rid <- randomForestSRC::var.select(object = fit, conservative = "high")
  return(rid$topvars)
}

#' Predict with RSF model
#'
#' @param fit RSF model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_rsf <- function(fit, newdata) {
  rs <- stats::predict(fit, newdata = newdata)$predicted
  return(round(as.numeric(rs), 10))
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
#' @noRd
train_enet <- function(est_dd, rid = NULL, alpha = 0.5, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x1, x2,
    family = "cox",
    alpha = alpha,
    nfolds = resolve_survival_cv_folds(
      est_dd$OS,
      requested_folds = 10L,
      context = "Enet Cox cross-validation"
    )
  )
  attr(fit, "iklsurvml_features") <- rid
  return(fit)
}

#' Predict with Enet model
#'
#' @param fit cv.glmnet model
#' @param newdata New data for prediction
#' @param rid Feature names to use
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_enet <- function(fit, newdata, rid) {
  return(as.numeric(stats::predict(
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
#' @noRd
train_lasso <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x1, x2,
    nfolds = resolve_survival_cv_folds(
      est_dd$OS,
      requested_folds = 10L,
      context = "Lasso Cox cross-validation"
    ),
    family = "cox",
    alpha = 1
  )
  attr(fit, "iklsurvml_features") <- rid
  return(fit)
}

#' Get Lasso selected variables
#'
#' @param fit cv.glmnet model
#' @return Character vector of selected variable names
#' @keywords internal
#' @noRd
get_lasso_selected_vars <- function(fit) {
  my_coefs <- stats::coef(fit, s = "lambda.min")
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
#' @noRd
predict_lasso <- function(fit, newdata, rid) {
  return(as.numeric(stats::predict(
    fit,
    type = "link",
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
#' @noRd
train_ridge <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  x1 <- as.matrix(est_dd[, rid])
  x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
  set.seed(seed)
  # Ridge regression: alpha = 0
  fit <- glmnet::glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
  cv_fit <- glmnet::cv.glmnet(
    x1, x2,
    nfolds = resolve_survival_cv_folds(
      est_dd$OS,
      requested_folds = 10L,
      context = "Ridge Cox cross-validation"
    ),
    family = "cox",
    alpha = 0
  )
  attr(fit, "iklsurvml_features") <- rid
  out <- list(fit = fit, cv.fit = cv_fit)
  attr(out, "iklsurvml_features") <- rid
  return(out)
}

#' Predict with Ridge model
#'
#' @param fit cv.glmnet model
#' @param newdata New data for prediction
#' @param rid Feature names to use
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
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
  return(as.numeric(stats::predict(
    glmnet_fit,
    type = "link",
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
#' @noRd
train_stepcox <- function(est_dd, direction = "both") {
  feature_names <- setdiff(colnames(est_dd), c("OS.time", "OS"))
  full_formula <- stats::reformulate(feature_names, response = "survival::Surv(OS.time, OS)")
  null_formula <- survival::Surv(OS.time, OS) ~ 1
  if (identical(direction, "forward")) {
    fit <- stats::step(
      survival::coxph(null_formula, data = est_dd),
      scope = list(
        lower = null_formula,
        upper = full_formula
      ),
      direction = direction,
      trace = 0
    )
    attr(fit, "iklsurvml_features") <- names(stats::coef(fit))
    return(fit)
  }

  fit <- stats::step(
    survival::coxph(full_formula, data = est_dd),
    direction = direction,
    trace = 0
  )
  attr(fit, "iklsurvml_features") <- names(stats::coef(fit))
  return(fit)
}

#' Get StepCox selected variables
#'
#' @param fit coxph model from stepcox
#' @return Character vector of selected variable names
#' @keywords internal
#' @noRd
get_stepcox_selected_vars <- function(fit) {
  return(names(stats::coef(fit)))
}

#' Predict with StepCox model
#'
#' @param fit coxph model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_stepcox <- function(fit, newdata) {
  return(stats::predict(fit, type = "risk", newdata = newdata))
}

# ---- CoxBoost ----

#' Train CoxBoost model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @return Trained CoxBoost model
#' @keywords internal
#' @noRd
train_coxboost <- function(est_dd, seed = 5201314) {
  set.seed(seed)
  features <- colnames(est_dd)[-c(1, 2)]
  pen <- CoxBoost::optimCoxBoostPenalty(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    trace = FALSE,
    start.penalty = 500,
    parallel = FALSE  # 禁用并行以确保可重复性
  )

  set.seed(seed)
  cv_res <- CoxBoost::cv.CoxBoost(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    maxstepno = 500,
    K = resolve_survival_cv_folds(
      est_dd$OS,
      requested_folds = 10L,
      context = "CoxBoost cross-validation"
    ),
    type = "verweij",
    penalty = pen$penalty
  )

  set.seed(seed)
  fit <- CoxBoost::CoxBoost(
    est_dd[, "OS.time"],
    est_dd[, "OS"],
    as.matrix(est_dd[, -c(1, 2)]),
    stepno = cv_res$optimal.step,
    penalty = pen$penalty
  )
  attr(fit, "iklsurvml_features") <- features
  return(fit)
}

#' Get CoxBoost selected variables
#'
#' @param fit CoxBoost model
#' @param tol Coefficient tolerance for non-zero selection
#' @return Character vector of selected variable names
#' @keywords internal
#' @noRd
get_coxboost_selected_vars <- function(fit, tol = .Machine$double.eps^0.5) {
  coefs <- stats::coef(fit)
  selected <- names(coefs)[!is.na(coefs) & abs(coefs) > tol]
  return(selected)
}

#' Predict with CoxBoost model
#'
#' @param fit CoxBoost model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_coxboost <- function(fit, newdata) {
  predictor_data <- newdata
  if (!is.null(fit[["xnames"]])) {
    missing <- setdiff(fit[["xnames"]], colnames(newdata))
    if (length(missing) > 0L) {
      stop(paste0("CoxBoost prediction data is missing columns: ",
                  paste(missing, collapse = ", ")), call. = FALSE)
    }
    predictor_data <- newdata[, fit[["xnames"]], drop = FALSE]
  } else if (all(c("OS.time", "OS") %in% colnames(newdata))) {
    predictor_data <- newdata[, setdiff(colnames(newdata), c("OS.time", "OS")), drop = FALSE]
  }
  return(as.numeric(stats::predict(
    fit,
    newdata = predictor_data,
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
#' @noRd
train_plsrcox <- function(est_dd, rid = NULL, seed = 5201314) {
  if (is.null(rid)) {
    rid <- colnames(est_dd)[-c(1, 2)]
  }
  prep <- prepare_plsrcox_features(est_dd, rid)
  selected_nt <- select_plsrcox_component_count(
    prep = prep,
    time = est_dd$OS.time,
    status = est_dd$OS,
    seed = seed
  )
  fit <- fit_plsrcox_with_fallback(
    prep = prep,
    time = est_dd$OS.time,
    status = est_dd$OS,
    selected_nt = selected_nt
  )
  attr(fit, "iklsurvml_features") <- rid
  attr(fit, "iklsurvml_center") <- prep$center
  attr(fit, "iklsurvml_scale") <- prep$scale
  attr(fit, "iklsurvml_selected_nt") <- selected_nt
  return(fit)
}

#' Select a stable plsRcox component count
#'
#' plsRcox cross-validation can be fold-partition sensitive on small or
#' near-separable survival data. Try the historical 5-fold setting first,
#' then fall back to alternate fold counts before using one safe component.
#' @keywords internal
#' @noRd
select_plsrcox_component_count <- function(prep, time, status, seed,
                                           max_nt = min(10, ncol(prep$x)),
                                           cv_fun = plsRcox::cv.plsRcox) {
  max_nt <- max(1L, as.integer(max_nt))
  n_obs <- length(time)
  nfold_candidates <- unique(pmin(c(5L, 10L, 3L), n_obs))
  nfold_candidates <- nfold_candidates[nfold_candidates >= 2L]
  if (length(nfold_candidates) == 0L) {
    return(1L)
  }

  errors <- character()
  for (nfold in nfold_candidates) {
    set.seed(seed)
    cv_res <- tryCatch(
      suppressWarnings(cv_fun(
        list(
          x = prep$x,
          time = time,
          status = status
        ),
        nt = max_nt,
        nfold = nfold,
        plot.it = FALSE,
        verbose = FALSE
      )),
      error = function(e) e
    )
    if (!inherits(cv_res, "error")) {
      selected_nt <- normalize_plsrcox_component_count(cv_res, max_nt)
      if (!is.na(selected_nt)) {
        return(selected_nt)
      }
      errors <- c(errors, sprintf("%d-fold returned no valid component count", nfold))
    } else {
      errors <- c(errors, sprintf("%d-fold: %s", nfold, conditionMessage(cv_res)))
    }
  }

  warning(
    "plsRcox cross-validation failed; falling back to 1 component",
    if (length(errors) > 0L) paste0(" (", paste(unique(errors), collapse = "; "), ")") else "",
    call. = FALSE
  )
  1L
}

#' Normalize plsRcox cross-validation output into a component count
#' @keywords internal
#' @noRd
normalize_plsrcox_component_count <- function(cv_res, max_nt) {
  candidates <- c(
    if (is.list(cv_res) && length(cv_res) >= 5L) cv_res[[5L]] else NA,
    if (is.list(cv_res) && !is.null(cv_res$nt)) cv_res$nt else NA,
    if (is.list(cv_res) && !is.null(cv_res$ncomp)) cv_res$ncomp else NA
  )
  candidates <- suppressWarnings(as.integer(candidates))
  candidates <- candidates[is.finite(candidates) & !is.na(candidates) & candidates >= 1L]
  if (length(candidates) == 0L) {
    return(NA_integer_)
  }
  min(candidates[[1L]], max_nt)
}

#' Fit plsRcox and reduce components if the final Cox fit is unstable
#' @keywords internal
#' @noRd
fit_plsrcox_with_fallback <- function(prep, time, status, selected_nt) {
  max_nt <- max(1L, min(as.integer(selected_nt), ncol(prep$x)))
  component_candidates <- unique(c(max_nt, rev(seq_len(max_nt))))
  errors <- character()

  for (nt in component_candidates) {
    fit <- tryCatch(
      suppressWarnings(plsRcox::plsRcox(
        prep$x,
        time = time,
        event = status,
        nt = nt,
        scaleX = FALSE,
        verbose = FALSE
      )),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      attr(fit, "iklsurvml_fit_nt") <- nt
      return(fit)
    }
    errors <- c(errors, sprintf("%d component(s): %s", nt, conditionMessage(fit)))
  }

  stop(
    "plsRcox training failed after component fallback: ",
    paste(unique(errors), collapse = "; "),
    call. = FALSE
  )
}

#' Standardize plsRcox feature matrices for numeric stability
#' @keywords internal
#' @noRd
prepare_plsrcox_features <- function(data, rid, center = NULL, scale = NULL) {
  x <- as.data.frame(data[, rid, drop = FALSE])
  x[] <- lapply(x, as.numeric)
  mat <- as.matrix(x)

  if (is.null(center)) {
    center <- colMeans(mat, na.rm = TRUE)
  }
  if (is.null(scale)) {
    scale <- apply(mat, 2, stats::sd, na.rm = TRUE)
  }
  scale[is.na(scale) | scale == 0] <- 1

  scaled <- sweep(mat, 2, center, FUN = "-")
  scaled <- sweep(scaled, 2, scale, FUN = "/")
  scaled <- as.data.frame(scaled, check.names = FALSE)
  colnames(scaled) <- rid

  list(x = scaled, center = center, scale = scale)
}

#' Predict with plsRcox model
#'
#' @param fit plsRcox model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_plsrcox <- function(fit, newdata) {
  rid <- attr(fit, "iklsurvml_features")
  if (is.null(rid)) {
    rid <- colnames(fit[["dataX"]])
  }
  center <- attr(fit, "iklsurvml_center")
  scale <- attr(fit, "iklsurvml_scale")
  if (!is.null(center) && !is.null(scale)) {
    newdata <- prepare_plsrcox_features(newdata, rid, center, scale)$x
  } else {
    newdata <- newdata[, rid, drop = FALSE]
  }
  return(as.numeric(stats::predict(
    fit,
    type = "lp",
    newdata = newdata
  )))
}

# ---- SuperPC ----

#' Train SuperPC model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @return List with fit and cv.fit
#' @keywords internal
#' @noRd
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

  # Retry on transient CV errors
  cv_fit <- NULL
  max_attempts <- 3
  for (i in seq_len(max_attempts)) {
    tryCatch({
      set.seed(seed)  # 确保交叉验证的可重复性
      cv_fit <- superpc::superpc.cv(
        fit, data,
        n.threshold = 20,
        n.fold = resolve_survival_cv_folds(
          est_dd$OS,
          requested_folds = 10L,
          context = "SuperPC cross-validation"
        ),
        n.components = 3,
        min.features = 2,
        max.features = nrow(data$x),
        compute.fullcv = TRUE,
        compute.preval = TRUE
      )
      break
    }, error = function(e) {
      warning(paste0("SuperPC CV attempt ", i, " failed: ", conditionMessage(e)))
      Sys.sleep(1)
    })
  }

  if (is.null(cv_fit)) {
    stop("SuperPC cross-validation failed after all retry attempts")
  }

  return(make_superpc_model(fit = fit, cv_fit = cv_fit, train_data = est_dd))
}

#' Build a SuperPC model object with its training feature subset
#'
#' @keywords internal
#' @noRd
make_superpc_model <- function(fit, cv_fit, train_data) {
  features <- colnames(train_data)[-c(1, 2)]
  train_data <- as.data.frame(train_data[, c("OS.time", "OS", features), drop = FALSE])
  prediction_params <- select_superpc_prediction_params(cv_fit)
  model <- structure(
    list(
      fit = fit,
      cv_fit = cv_fit,
      train_data = train_data,
      features = features,
      prediction_params = prediction_params
    ),
    class = c("ikl_superpc_model", "list")
  )
  attr(model, "iklsurvml_features") <- features
  model
}

#' Select SuperPC prediction threshold and component count from CV scores
#' @keywords internal
#' @noRd
select_superpc_prediction_params <- function(cv_fit) {
  thresholds <- cv_fit$thresholds
  scor <- cv_fit[["scor"]]
  if (is.null(thresholds) || length(thresholds) == 0L || is.null(scor)) {
    fallback_threshold <- if (!is.null(thresholds) && length(thresholds) > 0L) thresholds[[1L]] else 0
    return(list(threshold = fallback_threshold, n.components = 1L))
  }

  if (is.matrix(scor) || length(dim(scor)) == 2L) {
    best <- which(scor == max(scor, na.rm = TRUE), arr.ind = TRUE)[1, , drop = FALSE]
    component <- as.integer(best[1, 1])
    threshold_index <- as.integer(best[1, 2])
  } else {
    threshold_index <- which.max(scor)
    component <- 1L
  }

  threshold_index <- max(1L, min(threshold_index, length(thresholds)))
  component <- max(1L, component)
  list(
    threshold = thresholds[[threshold_index]],
    n.components = component
  )
}

#' Extract SuperPC model components from current or legacy objects
#'
#' @keywords internal
#' @noRd
extract_superpc_model <- function(superpc_obj) {
  if (!is.list(superpc_obj)) {
    stop("SuperPC model must be a list-like object", call. = FALSE)
  }

  fit <- if (!is.null(superpc_obj$fit)) superpc_obj$fit else superpc_obj[[1]]
  cv_fit <- if (!is.null(superpc_obj$cv_fit)) {
    superpc_obj$cv_fit
  } else if (!is.null(superpc_obj$cv.fit)) {
    superpc_obj$cv.fit
  } else {
    superpc_obj[[2]]
  }
  train_data <- if (!is.null(superpc_obj$train_data)) {
    superpc_obj$train_data
  } else if (!is.null(superpc_obj$training_data)) {
    superpc_obj$training_data
  } else {
    NULL
  }
  features <- superpc_obj$features
  prediction_params <- superpc_obj$prediction_params

  if (is.null(features) && !is.null(train_data)) {
    features <- colnames(train_data)[-c(1, 2)]
  }
  if (is.null(fit) || is.null(cv_fit)) {
    stop("SuperPC model is missing fit or cv_fit", call. = FALSE)
  }

  if (is.null(prediction_params)) {
    prediction_params <- select_superpc_prediction_params(cv_fit)
  }

  list(
    fit = fit,
    cv_fit = cv_fit,
    train_data = train_data,
    features = features,
    prediction_params = prediction_params
  )
}

#' Predict with SuperPC model
#'
#' @param fit SuperPC model
#' @param cv_fit Cross-validation result
#' @param train_data Training data used for model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
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

  prediction_params <- select_superpc_prediction_params(cv_fit)
  ff <- superpc::superpc.predict(
    fit, data, test,
    threshold = prediction_params$threshold,
    n.components = prediction_params$n.components
  )
  v_pred <- ff$v.pred
  n_obs <- nrow(newdata)
  component <- as.integer(prediction_params$n.components)

  if (is.matrix(v_pred) || is.data.frame(v_pred)) {
    v_pred <- as.matrix(v_pred)
    if (nrow(v_pred) == n_obs) {
      component <- max(1L, min(component, ncol(v_pred)))
      risk <- v_pred[, component]
    } else if (ncol(v_pred) == n_obs) {
      component <- max(1L, min(component, nrow(v_pred)))
      risk <- v_pred[component, ]
    } else {
      stop(paste0(
        "SuperPC prediction returned a ", paste(dim(v_pred), collapse = "x"),
        " matrix for ", n_obs, " rows; cannot derive one risk score per sample"
      ), call. = FALSE)
    }
  } else {
    risk <- as.numeric(v_pred)
  }

  if (length(risk) != n_obs) {
    stop(paste0(
      "SuperPC prediction returned ", length(risk),
      " risk scores for ", n_obs, " rows"
    ), call. = FALSE)
  }
  return(round(as.numeric(risk), 10))
}

#' Predict with a wrapped SuperPC model while preserving its feature subset
#'
#' @keywords internal
#' @noRd
predict_superpc_model <- function(superpc_obj, fallback_train_data, newdata) {
  parts <- extract_superpc_model(superpc_obj)
  train_data <- parts$train_data
  if (is.null(train_data)) {
    train_data <- fallback_train_data
  }
  features <- parts$features
  if (is.null(features)) {
    features <- colnames(train_data)[-c(1, 2)]
  }

  required <- c("OS.time", "OS", features)
  missing_train <- setdiff(required, colnames(train_data))
  missing_new <- setdiff(required, colnames(newdata))
  if (length(missing_train) > 0) {
    stop(paste0("SuperPC training data is missing columns: ",
                paste(missing_train, collapse = ", ")), call. = FALSE)
  }
  if (length(missing_new) > 0) {
    stop(paste0("SuperPC prediction data is missing columns: ",
                paste(missing_new, collapse = ", ")), call. = FALSE)
  }

  train_subset <- as.data.frame(train_data[, required, drop = FALSE])
  new_subset <- as.data.frame(newdata[, required, drop = FALSE])
  predict_superpc(parts$fit, parts$cv_fit, train_subset, new_subset)
}

# ---- GBM ----

#' Train GBM model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @param cores_for_parallel Number of cores for parallel processing
#' @return List with fit and best iteration
#' @keywords internal
#' @noRd
train_gbm <- function(est_dd, seed = 5201314, cores_for_parallel = 6) {
  # gbm's internal CV prediction path for distribution="coxph" has been
  # observed to segfault on small survival datasets in R 4.5/gbm builds.
  # Avoid that native-code path entirely and use OOB/train-error selection
  # from a single-process fit instead.
  n_trees <- 10000L
  bag_fraction <- 0.8
  minobs <- as.integer(max(
    1L,
    min(10L, floor((nrow(est_dd) * bag_fraction - 2) / 2))
  ))
  set.seed(seed)
  fit <- gbm::gbm(
    formula = survival::Surv(OS.time, OS) ~ .,
    data = est_dd,
    distribution = "coxph",
    n.trees = n_trees,
    interaction.depth = 3,
    n.minobsinnode = minobs,
    shrinkage = 0.001,
    bag.fraction = bag_fraction,
    cv.folds = 0,
    n.cores = 1
  )

  # Find optimal number of trees without invoking gbmCrossValPredictions.
  selection_method <- "OOB"
  best <- tryCatch(
    suppressMessages(suppressWarnings(gbm::gbm.perf(fit, method = "OOB", plot.it = FALSE))),
    error = function(e) NA_integer_
  )
  if (length(best) != 1 || is.na(best) || best < 1) {
    selection_method <- "train.error"
    best <- which.min(fit$train.error)
  }
  if (length(best) != 1 || is.na(best) || best < 1) {
    selection_method <- "last"
    best <- n_trees
  }
  best <- as.integer(best)

  set.seed(seed)
  fit <- gbm::gbm(
    formula = survival::Surv(OS.time, OS) ~ .,
    data = est_dd,
    distribution = "coxph",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = minobs,
    shrinkage = 0.001,
    bag.fraction = bag_fraction,
    cv.folds = 0,
    n.cores = 1
  )

  features <- colnames(est_dd)[-c(1, 2)]
  attr(fit, "iklsurvml_features") <- features
  out <- list(fit = fit, best = best, selection_method = selection_method)
  attr(out, "iklsurvml_features") <- features
  attr(out, "iklsurvml_gbm_selection_method") <- selection_method
  return(out)
}

#' Predict with GBM model
#'
#' @param fit GBM model
#' @param best Best iteration number
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_gbm <- function(fit, best, newdata) {
  return(as.numeric(stats::predict(fit, newdata, n.trees = best, type = "link")))
}

# ---- survivalsvm ----

#' Train survivalsvm model
#'
#' @param est_dd Training data with OS.time, OS and features
#' @param seed Random seed
#' @return Trained survivalsvm model
#' @keywords internal
#' @noRd
train_survivalsvm <- function(est_dd, seed = 5201314) {
  attempts <- list(
    list(gamma.mu = 1, opt.meth = "quadprog"),
    list(gamma.mu = 1, opt.meth = "ipop"),
    list(gamma.mu = 10, opt.meth = "quadprog"),
    list(gamma.mu = 0.1, opt.meth = "quadprog")
  )
  errors <- character()

  for (attempt in attempts) {
    set.seed(seed)
    fit <- tryCatch(
      survivalsvm::survivalsvm(
        survival::Surv(OS.time, OS) ~ .,
        data = est_dd,
        gamma.mu = attempt$gamma.mu,
        opt.meth = attempt$opt.meth
      ),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      attr(fit, "iklsurvml_gamma.mu") <- attempt$gamma.mu
      attr(fit, "iklsurvml_opt.meth") <- attempt$opt.meth
      attr(fit, "iklsurvml_features") <- colnames(est_dd)[-c(1, 2)]
      return(fit)
    }
    errors <- c(errors, paste0(
      "gamma.mu=", attempt$gamma.mu,
      ", opt.meth=", attempt$opt.meth,
      ": ", conditionMessage(fit)
    ))
  }

  stop(
    "survival-SVM training failed after solver/penalty fallbacks: ",
    paste(errors, collapse = "; "),
    call. = FALSE
  )
}

#' Predict with survivalsvm model
#'
#' @param fit survivalsvm model
#' @param newdata New data for prediction
#' @return Predicted risk scores
#' @keywords internal
#' @noRd
predict_survivalsvm <- function(fit, newdata) {
  # survivalsvm's default regression output is survival-time oriented on the
  # package's benchmark data: larger predictions indicate better prognosis.
  # iklSurvML risk scores use the opposite convention everywhere else
  # (larger RS = higher risk), so invert here and route all downstream
  # recalculation through this wrapper.
  return(-as.numeric(stats::predict(fit, newdata)$predicted))
}

# ---- Utility functions ----

#' Calculate C-index for predictions
#'
#' @param rs Risk scores
#' @param data Data with OS.time and OS
#' @return C-index value
#' @keywords internal
#' @noRd
calculate_cindex <- function(rs, data) {
  cindex_data <- data.frame(
    OS.time = data$OS.time,
    OS = data$OS,
    rs = as.numeric(rs)
  )
  cindex_data <- cindex_data[stats::complete.cases(cindex_data), ]

  if (nrow(cindex_data) < 2 || length(unique(cindex_data$rs)) < 2) {
    return(NA_real_)
  }

  return(as.numeric(survival::concordance(
    survival::Surv(OS.time, OS) ~ rs,
    data = cindex_data,
    reverse = TRUE
  )$concordance))
}

#' Calculate risk scores for validation data
#'
#' @param val_dd_list List of validation data
#' @param rs_func Function to calculate risk scores
#' @return List of risk score data frames
#' @keywords internal
#' @noRd
calculate_risk_scores <- function(val_dd_list, rs_func) {
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = rs_func(x))
  })
  return(rs)
}
