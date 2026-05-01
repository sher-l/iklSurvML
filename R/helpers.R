# Helper functions for Mime package
#
# This file contains utility functions used across the package

#' Return ID column to risk score table
#'
#' @param rs.table.list A list of risk score tables
#' @param rawtableID A list of original tables with ID column
#' @return A list of risk score tables with ID column added
#' @keywords internal
return_id_to_rs <- function(rs.table.list, rawtableID) {
  for (i in names(rs.table.list)) {
    rs.table.list[[i]]$ID <- rawtableID[[i]]$ID
    # Reorder columns: ID first
    cols <- c("ID", setdiff(colnames(rs.table.list[[i]]), "ID"))
    rs.table.list[[i]] <- rs.table.list[[i]][, cols]
  }
  return(rs.table.list)
}

#' Normalize ML feature and column names
#'
#' Survival model development historically converted dashes to dots in gene
#' names.  Univariate Cox filtering also converted underscores, which meant
#' filtered feature names could no longer be found in the training data.  Keep
#' all ML entry points on the same reversible-enough convention.
#'
#' @param x Character vector of feature names
#' @return Normalized feature names
#' @keywords internal
normalize_ml_feature_names <- function(x) {
  gsub("[-_]", ".", x)
}

#' Normalize a data frame's column names for ML workflows
#' @keywords internal
normalize_ml_data_columns <- function(x, context = "data") {
  colnames(x) <- normalize_ml_feature_names(colnames(x))
  duplicated_cols <- colnames(x)[duplicated(colnames(x))]
  if (length(duplicated_cols) > 0) {
    stop(paste0(
      "Duplicate columns after feature-name normalization in ", context, ": ",
      paste(unique(duplicated_cols), collapse = ", ")
    ))
  }
  x
}

#' Safely coerce a vector to numeric for ML workflows
#'
#' Factors are first converted to character so factor levels are not silently
#' re-encoded as 1, 2, 3, ... .  Non-empty values that become NA are rejected
#' instead of being hidden by downstream mean imputation.
#' @keywords internal
safe_as_numeric <- function(x, column, label) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.logical(x) || is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }
  if (!is.character(x)) {
    stop(paste0(label, " column '", column, "' must be numeric or character-numeric; got class ",
                paste(class(x), collapse = "/")), call. = FALSE)
  }

  empty_or_missing <- is.na(x) | trimws(x) == ""
  out <- suppressWarnings(as.numeric(x))
  introduced_na <- is.na(out) & !empty_or_missing
  if (any(introduced_na)) {
    bad <- unique(x[introduced_na])
    bad <- bad[seq_len(min(length(bad), 5L))]
    stop(paste0(
      label, " column '", column, "' contains non-numeric values: ",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }
  out
}

#' Safely coerce selected data frame columns to numeric
#' @keywords internal
coerce_numeric_columns <- function(x, columns, label) {
  for (column in columns) {
    x[[column]] <- safe_as_numeric(x[[column]], column, label)
  }
  x
}

#' Validate that survival status uses the package-wide 0/1 event convention
#' @keywords internal
assert_binary_survival_status <- function(x, label = "data") {
  status_values <- sort(unique(stats::na.omit(x$OS)))
  if (length(status_values) == 0L) {
    stop(paste0(label, " has no non-missing OS status values"), call. = FALSE)
  }
  if (!all(status_values %in% c(0, 1))) {
    stop(paste0(
      label,
      " OS must be coded as 0 = censored/alive and 1 = event/dead; observed values: ",
      paste(status_values, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate sample IDs before model fitting or prediction
#' @keywords internal
assert_unique_ids <- function(x, label = "data", id_col = "ID") {
  if (!id_col %in% colnames(x)) {
    stop(paste0(label, " is missing required ID column: ", id_col), call. = FALSE)
  }
  ids <- as.character(x[[id_col]])
  if (any(is.na(ids) | !nzchar(trimws(ids)))) {
    stop(paste0(label, " contains missing or blank sample IDs"), call. = FALSE)
  }
  duplicated_ids <- unique(ids[duplicated(ids)])
  if (length(duplicated_ids) > 0L) {
    duplicated_ids <- duplicated_ids[seq_len(min(length(duplicated_ids), 5L))]
    stop(paste0(
      label, " contains duplicated sample IDs: ",
      paste(duplicated_ids, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve candidate features across training and validation cohorts
#'
#' The default strict mode prevents validation cohorts from silently shrinking
#' the training feature space.  The legacy intersection behavior remains
#' available through feature_alignment = "intersection".
#' @keywords internal
resolve_survival_common_features <- function(train_data,
                                             list_train_vali_Data,
                                             candidate_genes,
                                             feature_alignment = c("strict", "intersection")) {
  feature_alignment <- match.arg(feature_alignment)
  required_cols <- c("ID", "OS.time", "OS")
  candidate_genes <- normalize_ml_feature_names(candidate_genes)
  requested <- c(required_cols, candidate_genes)

  if (identical(feature_alignment, "strict")) {
    datasets <- c(list(`Training data` = train_data), list_train_vali_Data)
    missing_by_dataset <- vapply(names(datasets), function(nm) {
      paste(setdiff(requested, colnames(datasets[[nm]])), collapse = ", ")
    }, character(1))
    missing_by_dataset <- missing_by_dataset[nzchar(missing_by_dataset)]
    if (length(missing_by_dataset) > 0L) {
      stop(paste0(
        "feature_alignment='strict' requires ID, OS.time, OS, and all candidate_genes ",
        "in every training/validation dataset. Missing: ",
        paste(paste0(names(missing_by_dataset), " [", missing_by_dataset, "]"), collapse = "; ")
      ), call. = FALSE)
    }
    return(requested)
  }

  common_feature <- requested
  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }
  dropped <- setdiff(candidate_genes, common_feature[-seq_along(required_cols)])
  if (length(dropped) > 0L) {
    warning(paste0(
      "feature_alignment='intersection' dropped candidate features absent from at least one cohort: ",
      paste(dropped, collapse = ", ")
    ), call. = FALSE)
  }
  common_feature
}

#' Normalize survival ML algorithm input aliases
#' @keywords internal
normalize_survival_ml_name <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  aliases <- c(
    "rsf" = "RSF",
    "enet" = "Enet",
    "stepcox" = "StepCox",
    "coxboost" = "CoxBoost",
    "plsrcox" = "plsRcox",
    "superpc" = "SuperPC",
    "gbm" = "GBM",
    "survivalsvm" = "survivalsvm",
    "survival-svm" = "survivalsvm",
    "survival - svm" = "survivalsvm",
    "survival svm" = "survivalsvm",
    "ridge" = "Ridge",
    "lasso" = "Lasso"
  )
  key <- tolower(trimws(x))
  normalized <- unname(aliases[key])
  ifelse(is.na(normalized), x, normalized)
}

#' Display name for survival ML models
#' @keywords internal
display_survival_ml_name <- function(x) {
  x <- normalize_survival_ml_name(x)
  if (identical(x, "survivalsvm")) {
    "survival-SVM"
  } else {
    x
  }
}

#' Survival ML algorithm registry
#'
#' This is the single source of truth for the public survival ML names and
#' first-stage selector eligibility.  Training/prediction dispatch still uses
#' the historical wrappers, but validation and combination guards should read
#' from this registry instead of hard-coded, drifting vectors.
#' @keywords internal
survival_ml_registry <- function() {
  data.frame(
    name = c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox",
             "SuperPC", "GBM", "survivalsvm", "Ridge", "Lasso"),
    display = c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox",
                "SuperPC", "GBM", "survival-SVM", "Ridge", "Lasso"),
    first_stage_selector = c(TRUE, FALSE, TRUE, TRUE, FALSE,
                             FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

#' Valid survival ML algorithm names
#' @keywords internal
survival_ml_names <- function() {
  survival_ml_registry()$name
}

#' Valid first-stage survival ML selector names
#' @keywords internal
survival_ml_first_stage_names <- function() {
  registry <- survival_ml_registry()
  registry$name[registry$first_stage_selector]
}

#' Reject survival ML self-combinations
#' @keywords internal
assert_no_survival_self_combination <- function(double_ml1, double_ml2) {
  if (!is.null(double_ml1) && !is.null(double_ml2) &&
      length(double_ml1) == 1 && length(double_ml2) == 1 &&
      !is.na(double_ml1) && !is.na(double_ml2) &&
      identical(double_ml1, double_ml2)) {
    stop(paste0(
      "Self-combinations are not supported: ", double_ml1, " + ", double_ml2,
      ". Use a distinct first-stage feature selector and second-stage model."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Normalize the fixed all-mode model grid
#' @keywords internal
normalize_all_mode_model_grid <- function(model_grid = "117") {
  if (is.null(model_grid)) {
    return("117")
  }
  if (length(model_grid) != 1L || anyNA(model_grid)) {
    stop("model_grid only supports '117'", call. = FALSE)
  }

  model_grid <- as.character(model_grid)
  if (!identical(model_grid, "117")) {
    stop("model_grid only supports '117'", call. = FALSE)
  }
  model_grid
}

#' Return StepCox directions used as first-stage selectors in all mode
#' @keywords internal
all_mode_stepcox_selector_dirs <- function(model_grid = "117") {
  normalize_all_mode_model_grid(model_grid)
  c("both", "backward", "forward")
}

#' Return CoxBoost second-stage algorithms for the all-mode grid
#' @keywords internal
all_mode_coxboost_second_stage_algorithms <- function(model_grid = "117") {
  normalize_all_mode_model_grid(model_grid)
  c("Enet", "GBM", "Lasso", "plsRcox", "Ridge", "StepCox", "SuperPC", "survivalsvm")
}

#' Return Lasso second-stage algorithms for the all-mode grid
#' @keywords internal
all_mode_lasso_second_stage_algorithms <- function(model_grid = "117") {
  normalize_all_mode_model_grid(model_grid)
  c("CoxBoost", "GBM", "plsRcox", "RSF", "StepCox", "SuperPC", "survivalsvm")
}

#' Return the Enet alpha values expanded by all-mode grids
#' @keywords internal
all_mode_alpha_values <- function() {
  seq(0.1, 0.9, by = 0.1)
}

#' Format one survival ML stage for model names
#' @keywords internal
format_survival_model_stage <- function(algorithm, parameter = NA_character_) {
  algorithm <- normalize_survival_ml_name(algorithm)
  if (identical(algorithm, "Enet")) {
    return(paste0("Enet[\u03b1=", parameter, "]"))
  }
  if (identical(algorithm, "StepCox")) {
    return(paste0("StepCox[", parameter, "]"))
  }
  display_survival_ml_name(algorithm)
}

#' Declarative fixed all-mode survival ML task table
#'
#' This table is the source of truth for the advertised 117-model grid.  Runner
#' implementations may still choose sequential, cached, or forked execution, but
#' completion checks compare their materialized model names against this table so
#' duplicated runner code cannot silently drift.
#' @keywords internal
survival_all_mode_task_table <- function(model_grid = "117") {
  normalize_all_mode_model_grid(model_grid)
  tasks <- list()

  add_task <- function(selector = NA_character_,
                       selector_param = NA_character_,
                       learner,
                       learner_param = NA_character_,
                       phase) {
    learner_name <- format_survival_model_stage(learner, learner_param)
    model_name <- learner_name
    if (!is.na(selector)) {
      model_name <- paste(
        format_survival_model_stage(selector, selector_param),
        learner_name,
        sep = " + "
      )
    }
    tasks[[length(tasks) + 1L]] <<- data.frame(
      phase = phase,
      selector = selector,
      selector_param = selector_param,
      learner = normalize_survival_ml_name(learner),
      learner_param = learner_param,
      model_name = model_name,
      stringsAsFactors = FALSE
    )
  }

  add_task(learner = "RSF", phase = "single")
  for (alpha in all_mode_alpha_values()) {
    add_task(learner = "Enet", learner_param = as.character(alpha), phase = "single")
  }
  for (direction in c("both", "backward", "forward")) {
    add_task(learner = "StepCox", learner_param = direction, phase = "single")
  }
  for (learner in c("CoxBoost", "plsRcox", "SuperPC", "GBM", "survivalsvm", "Ridge", "Lasso")) {
    add_task(learner = learner, phase = "single")
  }

  add_task(selector = "RSF", learner = "CoxBoost", phase = "RSF")
  for (alpha in all_mode_alpha_values()) {
    add_task(selector = "RSF", learner = "Enet",
             learner_param = as.character(alpha), phase = "RSF")
  }
  for (learner in c("GBM", "Lasso", "plsRcox", "Ridge", "SuperPC", "survivalsvm")) {
    add_task(selector = "RSF", learner = learner, phase = "RSF")
  }
  for (direction in c("both", "backward", "forward")) {
    add_task(selector = "RSF", learner = "StepCox",
             learner_param = direction, phase = "RSF")
  }

  for (selector_direction in all_mode_stepcox_selector_dirs(model_grid)) {
    add_task(selector = "StepCox", selector_param = selector_direction,
             learner = "CoxBoost", phase = "StepCox")
    for (alpha in all_mode_alpha_values()) {
      add_task(selector = "StepCox", selector_param = selector_direction,
               learner = "Enet", learner_param = as.character(alpha),
               phase = "StepCox")
    }
    for (learner in c("GBM", "Lasso", "plsRcox", "Ridge", "RSF", "SuperPC", "survivalsvm")) {
      add_task(selector = "StepCox", selector_param = selector_direction,
               learner = learner, phase = "StepCox")
    }
  }

  for (alpha in all_mode_alpha_values()) {
    add_task(selector = "CoxBoost", learner = "Enet",
             learner_param = as.character(alpha), phase = "CoxBoost")
  }
  for (learner in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add_task(selector = "CoxBoost", learner = learner, phase = "CoxBoost")
  }
  for (direction in c("both", "backward", "forward")) {
    add_task(selector = "CoxBoost", learner = "StepCox",
             learner_param = direction, phase = "CoxBoost")
  }
  for (learner in c("SuperPC", "survivalsvm")) {
    add_task(selector = "CoxBoost", learner = learner, phase = "CoxBoost")
  }

  for (learner in c("CoxBoost", "GBM", "plsRcox", "RSF")) {
    add_task(selector = "Lasso", learner = learner, phase = "Lasso")
  }
  for (direction in c("both", "backward", "forward")) {
    add_task(selector = "Lasso", learner = "StepCox",
             learner_param = direction, phase = "Lasso")
  }
  for (learner in c("SuperPC", "survivalsvm")) {
    add_task(selector = "Lasso", learner = learner, phase = "Lasso")
  }

  do.call(rbind, tasks)
}

#' Return the advertised all-mode model count for the fixed grid
#' @keywords internal
all_mode_model_grid_size <- function(model_grid = "117") {
  nrow(survival_all_mode_task_table(model_grid))
}

#' Return all expected all-mode model names in deterministic order
#' @keywords internal
all_mode_model_names <- function(model_grid = "117") {
  survival_all_mode_task_table(model_grid)$model_name
}

#' Warn when an all-mode run cannot materialize the full headline model grid
#'
#' Some first-stage selectors can legitimately return fewer than two features on
#' weak or low-dimensional data. In that case the affected two-stage models
#' cannot be fit, but already materialized models are still useful and should be
#' returned instead of being discarded by a final count assertion.
warn_if_all_mode_incomplete <- function(actual, expected = all_mode_model_grid_size(), context = "All-mode") {
  actual <- as.integer(actual)
  expected <- as.integer(expected)
  if (!identical(actual, expected)) {
    warning(paste0(
      context, " produced ", actual, " models; expected ", expected,
      ". Returning the successfully fitted models. Missing models usually mean ",
      "a first-stage feature selector returned fewer than 2 variables on this dataset."
    ), call. = FALSE)
  }
  invisible(actual)
}

#' Assert all-mode produced the complete advertised survival model grid
#'
#' Returning a partial all-mode result can make downstream model comparisons look
#' complete when they are not.  Keep an explicit opt-in escape hatch for weak
#' datasets whose first-stage selectors legitimately cannot materialize every
#' two-stage model.
#' @keywords internal
assert_complete_all_mode_result <- function(result,
                                            expected = all_mode_model_grid_size(),
                                            allow_partial = FALSE,
                                            context = "All-mode") {
  if (!is.logical(allow_partial) || length(allow_partial) != 1L ||
      is.na(allow_partial)) {
    stop("allow_partial must be TRUE or FALSE", call. = FALSE)
  }
  if (isTRUE(allow_partial)) {
    return(invisible(result))
  }

  actual <- if (!is.null(result$ml.res)) length(result$ml.res) else 0L
  model_errors <- if (!is.null(result$Model.errors)) result$Model.errors else character()
  model_skips <- if (!is.null(result$Model.skips)) result$Model.skips else character()

  if (!identical(as.integer(actual), as.integer(expected)) ||
      length(model_errors) > 0L ||
      length(model_skips) > 0L) {
    details <- c(
      paste0("built ", actual, " of ", expected, " expected models"),
      if (length(model_errors) > 0L) {
        paste0("errors: ", paste(model_errors, collapse = "; "))
      },
      if (length(model_skips) > 0L) {
        paste0("skips: ", paste(model_skips, collapse = "; "))
      }
    )
    stop(paste0(
      context,
      " did not complete the advertised model grid (",
      paste(details, collapse = " | "),
      "). Set allow_partial = TRUE to return partial results explicitly."
    ), call. = FALSE)
  }

  actual_names <- names(result$ml.res)
  expected_names <- all_mode_model_names()
  duplicate_names <- unique(actual_names[duplicated(actual_names)])
  missing_names <- setdiff(expected_names, actual_names)
  unexpected_names <- setdiff(actual_names, expected_names)
  if (length(duplicate_names) > 0L ||
      length(missing_names) > 0L ||
      length(unexpected_names) > 0L) {
    details <- c(
      if (length(duplicate_names) > 0L) {
        paste0("duplicate names: ", paste(duplicate_names, collapse = ", "))
      },
      if (length(missing_names) > 0L) {
        paste0("missing names: ", paste(missing_names, collapse = ", "))
      },
      if (length(unexpected_names) > 0L) {
        paste0("unexpected names: ", paste(unexpected_names, collapse = ", "))
      }
    )
    stop(paste0(
      context,
      " materialized a model set that does not match the declarative 117-grid (",
      paste(details, collapse = " | "),
      ")."
    ), call. = FALSE)
  }

  invisible(result)
}

#' Validate and canonicalize survival ML parameters
#' @keywords internal
validate_survival_ml_params <- function(mode,
                                        single_ml = NULL,
                                        double_ml1 = NULL,
                                        double_ml2 = NULL,
                                        alpha_for_enet = NULL,
                                        direction_for_stepcox = NULL) {
  valid_modes <- c("all", "single", "double")
  if (is.null(mode) || length(mode) != 1 || is.na(mode) || !mode %in% valid_modes) {
    stop(paste0("mode must be one of: ", paste(valid_modes, collapse = ", ")), call. = FALSE)
  }

  if (!is.null(alpha_for_enet) &&
      (length(alpha_for_enet) != 1 || is.na(alpha_for_enet) ||
       alpha_for_enet < 0 || alpha_for_enet > 1)) {
    stop("alpha_for_Enet must be between 0 and 1", call. = FALSE)
  }

  valid_directions <- c("both", "backward", "forward")
  if (is.null(direction_for_stepcox) ||
      length(direction_for_stepcox) != 1 ||
      is.na(direction_for_stepcox) ||
      !direction_for_stepcox %in% valid_directions) {
    stop(paste0(
      "direction_for_stepcox must be one of: ",
      paste(valid_directions, collapse = ", ")
    ), call. = FALSE)
  }

  single_ml <- normalize_survival_ml_name(single_ml)
  double_ml1 <- normalize_survival_ml_name(double_ml1)
  double_ml2 <- normalize_survival_ml_name(double_ml2)

  valid_single <- survival_ml_names()
  valid_first <- survival_ml_first_stage_names()

  if (identical(mode, "single") &&
      (is.null(single_ml) || length(single_ml) != 1 || is.na(single_ml) ||
       !single_ml %in% valid_single)) {
    stop(paste0("single_ml must be one of: ", paste(valid_single, collapse = ", ")), call. = FALSE)
  }

  if (identical(mode, "double")) {
    if (is.null(double_ml1) || length(double_ml1) != 1 || is.na(double_ml1) ||
        !double_ml1 %in% valid_first) {
      stop(paste0("double_ml1 must be one of: ", paste(valid_first, collapse = ", ")), call. = FALSE)
    }
    if (is.null(double_ml2) || length(double_ml2) != 1 || is.na(double_ml2) ||
        !double_ml2 %in% valid_single) {
      stop(paste0("double_ml2 must be one of: ", paste(valid_single, collapse = ", ")), call. = FALSE)
    }
    assert_no_survival_self_combination(double_ml1, double_ml2)
  }

  list(single_ml = single_ml, double_ml1 = double_ml1, double_ml2 = double_ml2)
}

#' Resolve requested model names from a stored survival ML result
#' @keywords internal
select_survival_model_names <- function(model_names,
                                        mode,
                                        single_ml = NULL,
                                        double_ml1 = NULL,
                                        double_ml2 = NULL,
                                        alpha_for_enet = NULL,
                                        direction_for_stepcox = NULL) {
  if (identical(mode, "all")) {
    return(model_names)
  }

  if (is.null(alpha_for_enet)) {
    alpha_for_enet <- 0.1
  }
  if (is.null(direction_for_stepcox)) {
    direction_for_stepcox <- "both"
  }

  algorithm_display <- function(algorithm) {
    algorithm <- normalize_survival_ml_name(algorithm)
    if (identical(algorithm, "Enet")) {
      paste0("Enet[\u03b1=", alpha_for_enet, "]")
    } else if (identical(algorithm, "StepCox")) {
      paste0("StepCox[", direction_for_stepcox, "]")
    } else {
      display_survival_ml_name(algorithm)
    }
  }

  if (identical(mode, "single")) {
    if (is.null(single_ml)) {
      stop("single_ml must be provided when mode='single'", call. = FALSE)
    }
    target <- algorithm_display(single_ml)
  } else if (identical(mode, "double")) {
    if (is.null(double_ml1) || is.null(double_ml2)) {
      stop("double_ml1 and double_ml2 must be provided when mode='double'", call. = FALSE)
    }
    assert_no_survival_self_combination(
      normalize_survival_ml_name(double_ml1),
      normalize_survival_ml_name(double_ml2)
    )
    target <- paste(algorithm_display(double_ml1), algorithm_display(double_ml2), sep = " + ")
  } else {
    stop("mode must be one of: all, single, double", call. = FALSE)
  }

  selected <- intersect(model_names, target)
  if (length(selected) == 0L) {
    stop(paste0(
      "Requested model '", target, "' was not found in res.by.ML.Dev.Prog.Sig$ml.res. ",
      "Available models: ", paste(model_names, collapse = ", ")
    ), call. = FALSE)
  }
  selected
}

#' Display progress bar for long operations
#'
#' @param index Current iteration index
#' @param totalN Total number of iterations
#' @param breakN Number of breaks in progress bar (default 20)
#' @keywords internal
display_progress <- function(index, totalN, breakN = 20) {
  if (index %% ceiling(totalN / breakN) == 0) {
    cat(paste(round(index * 100 / totalN), "% ", sep = ""))
  }
}

#' Preprocess data list for ML analysis
#'
#' @param list_train_vali_Data A list of training and validation data
#' @param common_feature Common features to select
#' @return Preprocessed data list
#' @keywords internal
preprocess_data_list <- function(list_train_vali_Data, common_feature, recipe = NULL) {
  data_names <- names(list_train_vali_Data)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(list_train_vali_Data))
  }
  list_train_vali_Data <- lapply(data_names, function(nm) {
    if (is.null(recipe)) {
      preprocess_ml_data_frame(list_train_vali_Data[[nm]], common_feature, paste0("Dataset '", nm, "'"))
    } else {
      apply_survival_preprocess_recipe(
        list_train_vali_Data[[nm]],
        recipe = recipe,
        label = paste0("Dataset '", nm, "'"),
        common_feature = common_feature
      )
    }
  })
  names(list_train_vali_Data) <- data_names
  list_train_vali_Data
}

#' Preprocess training data for ML analysis
#'
#' @param train_data Training data
#' @param common_feature Common features to select
#' @return Preprocessed training data
#' @keywords internal
preprocess_train_data <- function(train_data, common_feature, return_recipe = FALSE) {
  prepped <- fit_survival_preprocess_recipe(train_data, common_feature, "Training data")
  if (isTRUE(return_recipe)) {
    return(prepped)
  }
  prepped$data
}

#' Shared preprocessing for survival ML data frames
#' @keywords internal
preprocess_ml_data_frame <- function(x, common_feature, label) {
  fit_survival_preprocess_recipe(x, common_feature, label)$data
}

#' Fit a train-only survival preprocessing recipe
#' @keywords internal
fit_survival_preprocess_recipe <- function(x, common_feature, label = "Training data") {
  x <- normalize_ml_data_columns(x, label)
  common_feature <- normalize_ml_feature_names(common_feature)
  missing <- setdiff(common_feature, colnames(x))
  if (length(missing) > 0) {
    stop(paste0(label, " is missing required columns: ", paste(missing, collapse = ", ")), call. = FALSE)
  }

  x <- x[, common_feature, drop = FALSE]
  x <- coerce_numeric_columns(x, common_feature[-1], label)

  n_before <- nrow(x)
  x <- x[!is.na(x$OS.time) & !is.na(x$OS), , drop = FALSE]
  x <- x[x$OS.time > 0, , drop = FALSE]
  n_removed <- n_before - nrow(x)
  if (n_removed > 0) {
    warning(paste0(label, ": removed ", n_removed, " rows with invalid OS.time/OS"), call. = FALSE)
  }
  assert_binary_survival_status(x, label)
  assert_unique_ids(x, label)

  na_count <- sum(is.na(x[, -c(1:3), drop = FALSE]))
  feature_names <- common_feature[-c(1:3)]
  impute_values <- vapply(feature_names, function(feature) {
    col <- x[[feature]]
    if (all(is.na(col))) {
      stop(paste0(label, " has an all-NA feature column: ", feature), call. = FALSE)
    }
    mean(col, na.rm = TRUE)
  }, numeric(1))

  if (na_count > 0) {
    warning(paste0(label, ": ", na_count, " NA values imputed with column means"), call. = FALSE)
    for (feature in feature_names) {
      missing_idx <- is.na(x[[feature]])
      if (any(missing_idx)) {
        x[[feature]][missing_idx] <- impute_values[[feature]]
      }
    }
  }

  recipe <- structure(
    list(
      common_feature = common_feature,
      feature_names = feature_names,
      impute_values = impute_values
    ),
    class = c("ikl_survival_preprocess_recipe", "list")
  )

  list(data = x, recipe = recipe)
}

#' Apply a train-only survival preprocessing recipe
#' @keywords internal
apply_survival_preprocess_recipe <- function(x,
                                             recipe,
                                             label = "data",
                                             common_feature = NULL) {
  if (is.null(common_feature)) {
    common_feature <- recipe$common_feature
  }
  x <- normalize_ml_data_columns(x, label)
  common_feature <- normalize_ml_feature_names(common_feature)
  missing <- setdiff(common_feature, colnames(x))
  if (length(missing) > 0) {
    stop(paste0(label, " is missing required columns: ", paste(missing, collapse = ", ")), call. = FALSE)
  }

  x <- x[, common_feature, drop = FALSE]
  x <- coerce_numeric_columns(x, common_feature[-1], label)

  n_before <- nrow(x)
  x <- x[!is.na(x$OS.time) & !is.na(x$OS), , drop = FALSE]
  x <- x[x$OS.time > 0, , drop = FALSE]
  n_removed <- n_before - nrow(x)
  if (n_removed > 0) {
    warning(paste0(label, ": removed ", n_removed, " rows with invalid OS.time/OS"), call. = FALSE)
  }
  assert_binary_survival_status(x, label)
  assert_unique_ids(x, label)

  feature_names <- common_feature[-c(1:3)]
  missing_recipe_values <- setdiff(feature_names, names(recipe$impute_values))
  if (length(missing_recipe_values) > 0) {
    stop(paste0(
      "Training preprocessing recipe is missing imputation values for: ",
      paste(missing_recipe_values, collapse = ", ")
    ), call. = FALSE)
  }

  na_count <- sum(is.na(x[, feature_names, drop = FALSE]))
  if (na_count > 0) {
    warning(paste0(
      label, ": ", na_count,
      " NA values imputed with column means from training preprocessing recipe"
    ), call. = FALSE)
    for (feature in feature_names) {
      missing_idx <- is.na(x[[feature]])
      if (any(missing_idx)) {
        x[[feature]][missing_idx] <- recipe$impute_values[[feature]]
      }
    }
  }

  x
}

#' Validate and subset previous-signature input cohorts
#'
#' Missing genes in published prognostic signatures must not be silently
#' converted to expression value 0; absence from a cohort is different from
#' measured zero expression and changes the risk-score semantics.
#' @keywords internal
prepare_previous_signature_input_data <- function(list_input_data,
                                                  common_feature,
                                                  context = "previous signature") {
  data_names <- names(list_input_data)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(list_input_data))
  }

  missing_by_dataset <- vapply(data_names, function(nm) {
    missing <- setdiff(common_feature, colnames(list_input_data[[nm]]))
    paste(missing, collapse = ", ")
  }, character(1))
  missing_by_dataset <- missing_by_dataset[nzchar(missing_by_dataset)]

  if (length(missing_by_dataset) > 0L) {
    stop(paste0(
      context,
      " input data is missing required columns. Missing genes are not imputed as zero: ",
      paste(paste0(names(missing_by_dataset), " [", missing_by_dataset, "]"), collapse = "; ")
    ), call. = FALSE)
  }

  out <- lapply(data_names, function(nm) {
    list_input_data[[nm]][, common_feature, drop = FALSE]
  })
  names(out) <- data_names
  out
}

#' Strictly preprocess cohorts used by previous-signature survival helpers
#'
#' Previous-signature scores are fixed published formulas, so missing expression
#' values should not be cohort-wise mean-imputed or silently re-encoded.  Reject
#' unsafe inputs and require callers to provide complete numeric expression.
#' @keywords internal
preprocess_previous_signature_survival_data <- function(list_input_data,
                                                        common_feature,
                                                        context = "previous signature") {
  list_input_data <- prepare_previous_signature_input_data(
    list_input_data,
    common_feature,
    context = context
  )
  feature_names <- common_feature[-c(1:3)]

  data_names <- names(list_input_data)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(list_input_data))
  }
  out <- lapply(data_names, function(nm) {
    label <- paste0(context, " dataset '", nm, "'")
    x <- list_input_data[[nm]]
    x <- coerce_numeric_columns(x, c("OS.time", "OS", feature_names), label)

    n_before <- nrow(x)
    x <- x[!is.na(x$OS.time) & !is.na(x$OS), , drop = FALSE]
    x <- x[x$OS.time > 0, , drop = FALSE]
    n_removed <- n_before - nrow(x)
    if (n_removed > 0L) {
      warning(paste0(label, ": removed ", n_removed,
                     " rows with invalid OS.time/OS"), call. = FALSE)
    }
    assert_binary_survival_status(x, label)
    assert_unique_ids(x, label)

    na_count <- sum(is.na(x[, feature_names, drop = FALSE]))
    if (na_count > 0L) {
      stop(paste0(
        label,
        " contains ", na_count,
        " NA signature expression values; previous-signature scoring does not ",
        "impute missing expression values. Provide complete numeric data."
      ), call. = FALSE)
    }
    x
  })
  names(out) <- data_names
  out
}

#' Extract and score named immunotherapy signature genes
#'
#' This helper rejects absent genes instead of relying on positional slices after
#' `%in%` subsetting, which can silently shorten or reorder a signature matrix.
#' @keywords internal
calculate_signature_score_by_genes <- function(data,
                                               signature_genes,
                                               signature_name,
                                               dataset_name = "dataset",
                                               score_type = c("mean", "geometric", "single_gene"),
                                               single_gene = NULL) {
  score_type <- match.arg(score_type)
  label <- if (!is.null(dataset_name) && nzchar(dataset_name)) {
    paste0("Dataset '", dataset_name, "'")
  } else {
    "Dataset"
  }

  data <- normalize_ml_data_columns(data, label)
  signature_genes <- normalize_ml_feature_names(signature_genes)
  if (length(signature_genes) == 0L) {
    stop(paste0(signature_name, " has no signature genes"), call. = FALSE)
  }

  required <- unique(signature_genes)
  if (identical(score_type, "single_gene")) {
    if (is.null(single_gene) || length(single_gene) != 1L || is.na(single_gene)) {
      stop("single_gene must be provided when score_type='single_gene'", call. = FALSE)
    }
    single_gene <- normalize_ml_feature_names(single_gene)
    required <- unique(c(required, single_gene))
  }

  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop(paste0(
      signature_name, " in ", label,
      " is missing required genes: ",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }

  expr <- as.data.frame(data[, required, drop = FALSE], check.names = FALSE)
  expr <- coerce_numeric_columns(expr, required, paste0(label, " ", signature_name))

  if (identical(score_type, "mean")) {
    return(rowMeans(expr[, signature_genes, drop = FALSE]))
  }
  if (identical(score_type, "geometric")) {
    return(compositions::geometricmeanRow(expr[, signature_genes, drop = FALSE]))
  }

  as.numeric(expr[[single_gene]])
}

#' Calculate time-dependent ROC from fixed-direction risk scores
#' @keywords internal
calculate_survival_roc_from_risk <- function(x,
                                             AUC_time,
                                             auc_cal_method = "KM",
                                             risk_direction = "higher_is_worse") {
  auc_cal_method <- match.arg(auc_cal_method, c("KM", "NNE"))
  if (!risk_direction %in% c("higher_is_worse", "lower_is_worse")) {
    stop("risk_direction must be 'higher_is_worse' or 'lower_is_worse'", call. = FALSE)
  }

  x <- x[stats::complete.cases(x[, c("OS.time", "OS", "RS")]), , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("risk score table has no complete OS.time/OS/RS rows", call. = FALSE)
  }

  marker <- if (identical(risk_direction, "higher_is_worse")) x$RS else -x$RS
  roc_args <- list(
    Stime = x$OS.time,
    status = x$OS,
    marker = marker,
    predict.time = 365 * AUC_time,
    method = auc_cal_method
  )
  if (identical(auc_cal_method, "NNE")) {
    roc_args$span <- 0.25 * nrow(x)^(-0.20)
  }
  risk.survivalROC <- do.call(survivalROC::survivalROC, roc_args)

  HR <- NA_real_
  group <- ifelse(x$RS > stats::median(x$RS), "High", "Low")
  if (length(unique(group)) <= 1L) {
    group <- ifelse(x$RS > mean(x$RS), "High", "Low")
  }
  if (length(unique(group)) > 1L) {
    group <- factor(group, levels = c("Low", "High"))
    mySurv <- survival::Surv(x$OS.time, x$OS)
    data.survdiff <- survival::survdiff(mySurv ~ group)
    HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) /
      (data.survdiff$obs[1] / data.survdiff$exp[1])
  }

  roc_1 <- data.frame(
    TP = round(risk.survivalROC$TP, 3),
    FP = round(risk.survivalROC$FP, 3),
    AUC = risk.survivalROC$AUC,
    HR = HR,
    marker_direction = risk_direction
  )
  roc_1
}

#' Infer the final learner used by a survival model name
#' @keywords internal
survival_model_stage_algorithm <- function(stage_name) {
  if (grepl("^Enet\\[", stage_name)) {
    return("Enet")
  }
  if (grepl("^StepCox\\[", stage_name)) {
    return("StepCox")
  }
  normalize_survival_ml_name(stage_name)
}

#' Infer the final learner used by a survival model name
#' @keywords internal
survival_model_final_algorithm <- function(model_name) {
  final <- tail(strsplit(model_name, " \\+ ")[[1]], 1)
  survival_model_stage_algorithm(final)
}

#' Read iklSurvML metadata from raw or wrapped model objects
#' @keywords internal
get_iklsurvml_model_attr <- function(fit, attr_name) {
  value <- attr(fit, attr_name, exact = TRUE)
  if (!is.null(value)) {
    return(value)
  }
  if (is.list(fit) && !is.null(fit$fit)) {
    value <- attr(fit$fit, attr_name, exact = TRUE)
    if (!is.null(value)) {
      return(value)
    }
  }
  NULL
}

#' Extract the feature columns required to predict a stored survival model
#' @keywords internal
survival_model_features <- function(model_name, fit) {
  algorithm <- survival_model_final_algorithm(model_name)
  features <- get_iklsurvml_model_attr(fit, "iklsurvml_features")
  if (!is.null(features)) {
    return(features)
  }

  if (identical(algorithm, "RSF")) {
    return(fit[["xvar.names"]])
  }
  if (identical(algorithm, "survivalsvm")) {
    return(fit[["var.names"]])
  }
  if (identical(algorithm, "CoxBoost")) {
    return(fit[["xnames"]])
  }
  if (identical(algorithm, "Enet") || identical(algorithm, "Lasso")) {
    return(fit[["glmnet.fit"]][["beta"]]@Dimnames[[1]])
  }
  if (identical(algorithm, "plsRcox")) {
    features <- attr(fit, "iklsurvml_features")
    if (is.null(features)) {
      features <- colnames(fit[["dataX"]])
    }
    return(features)
  }
  if (identical(algorithm, "Ridge")) {
    ridge_fit <- if (is.list(fit) && "fit" %in% names(fit)) fit$fit else fit
    return(ridge_fit[["beta"]]@Dimnames[[1]])
  }
  if (identical(algorithm, "GBM")) {
    gbm_fit <- if (is.list(fit) && "fit" %in% names(fit)) fit$fit else fit[[1]]
    return(gbm_fit[["var.names"]])
  }
  if (identical(algorithm, "StepCox")) {
    return(names(stats::coef(fit)))
  }
  if (identical(algorithm, "SuperPC")) {
    return(extract_superpc_model(fit)$features)
  }

  stop(paste0("Unsupported survival model algorithm for: ", model_name), call. = FALSE)
}

#' Build a stable metadata record for a fitted survival model
#' @keywords internal
build_survival_model_info <- function(model_name, fit) {
  parts <- strsplit(model_name, " \\+ ")[[1]]
  learner <- survival_model_final_algorithm(model_name)
  selector <- if (length(parts) > 1L) survival_model_stage_algorithm(parts[[1L]]) else NULL
  features <- tryCatch(
    survival_model_features(model_name, fit),
    error = function(e) NULL
  )
  diagnostics <- list()
  if (identical(learner, "GBM") && is.list(fit)) {
    diagnostics$gbm_selection_method <- fit$selection_method
    diagnostics$gbm_best_trees <- fit$best
  }
  if (identical(learner, "plsRcox")) {
    diagnostics$plsrcox_selected_nt <- attr(fit, "iklsurvml_selected_nt", exact = TRUE)
    diagnostics$plsrcox_fit_nt <- attr(fit, "iklsurvml_fit_nt", exact = TRUE)
  }

  list(
    model_name = model_name,
    display_name = model_name,
    pipeline = list(selector = selector, learner = learner),
    features = features,
    risk_direction = "higher_is_worse",
    diagnostics = diagnostics
  )
}

#' Attach stable metadata records to survival ML result objects
#' @keywords internal
attach_survival_model_info <- function(result) {
  if (is.null(result) || is.null(result$ml.res) || length(result$ml.res) == 0L) {
    return(result)
  }
  result$Model.info <- stats::setNames(
    lapply(names(result$ml.res), function(model_name) {
      build_survival_model_info(model_name, result$ml.res[[model_name]])
    }),
    names(result$ml.res)
  )
  result
}

#' Subset survival prediction data to a model's required columns
#' @keywords internal
subset_survival_prediction_data <- function(x, features, model_name) {
  required <- c("OS.time", "OS", features)
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(paste0(
      model_name, " prediction data is missing columns: ",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  x[, required, drop = FALSE]
}

#' Recalculate risk scores for any stored survival model via one dispatch path
#' @keywords internal
calculate_survival_model_risk_scores <- function(model_name,
                                                 fit,
                                                 val_dd_list,
                                                 fallback_train_data,
                                                 model_info = NULL) {
  algorithm <- if (!is.null(model_info$pipeline$learner)) {
    normalize_survival_ml_name(model_info$pipeline$learner)
  } else {
    survival_model_final_algorithm(model_name)
  }
  features <- model_info$features
  if (is.null(features)) {
    features <- survival_model_features(model_name, fit)
  }
  if (is.null(features) && identical(algorithm, "SuperPC")) {
    features <- colnames(fallback_train_data)[-c(1, 2)]
  }

  if (identical(algorithm, "SuperPC")) {
    return(lapply(val_dd_list, function(x) {
      x_model <- subset_survival_prediction_data(x, features, model_name)
      cbind(x_model[, 1:2, drop = FALSE],
            RS = predict_superpc_model(fit, fallback_train_data, x_model))
    }))
  }

  val_model_list <- lapply(val_dd_list, subset_survival_prediction_data,
                           features = features, model_name = model_name)

  calculate_risk_scores(val_model_list, function(x) {
    if (identical(algorithm, "RSF")) {
      predict_rsf(fit, x)
    } else if (identical(algorithm, "survivalsvm")) {
      predict_survivalsvm(fit, x)
    } else if (identical(algorithm, "CoxBoost")) {
      predict_coxboost(fit, x)
    } else if (identical(algorithm, "Enet")) {
      predict_enet(fit, x, features)
    } else if (identical(algorithm, "Lasso")) {
      predict_lasso(fit, x, features)
    } else if (identical(algorithm, "plsRcox")) {
      predict_plsrcox(fit, x)
    } else if (identical(algorithm, "Ridge")) {
      predict_ridge(fit, x, features)
    } else if (identical(algorithm, "GBM")) {
      gbm_fit <- if (is.list(fit) && "fit" %in% names(fit)) fit$fit else fit[[1]]
      best <- if (is.list(fit) && "best" %in% names(fit)) fit$best else fit[[2]]
      predict_gbm(gbm_fit, best, x)
    } else if (identical(algorithm, "StepCox")) {
      predict_stepcox(fit, x)
    } else {
      stop(paste0("Unsupported survival model algorithm for: ", model_name),
           call. = FALSE)
    }
  })
}

#' Category ML optional engine dependency map
#' @keywords internal
category_method_packages <- function() {
  list(
    nb = "klaR",
    svmRadialWeights = "kernlab",
    rf = "randomForest",
    kknn = "kknn",
    adaboost = "fastAdaboost",
    LogitBoost = "caTools",
    cancerclass = c("cancerclass", "Biobase", "pROC")
  )
}

#' Category ML tune grids sized to the current feature count
#' @keywords internal
category_tune_grid <- function(n_features,
                               tune_profile = c("standard", "exhaustive")) {
  tune_profile <- match.arg(tune_profile)
  n_features <- max(1L, as.integer(n_features))
  rf_mtry <- unique(pmax(
    1L,
    pmin(n_features, round(seq(1, n_features, length.out = min(10L, n_features))))
  ))

  if (identical(tune_profile, "exhaustive")) {
    return(list(
      nb = expand.grid(fL = c(0, 0.5, 1, 1.5, 2.0),
                       usekernel = TRUE,
                       adjust = c(0.5, 0.75, 1, 1.25, 1.5)),
      svmRadialWeights = expand.grid(sigma = c(0.0005, 0.001, 0.005, 0.01, 0.05),
                                     C = c(1, 3, 5, 10, 20),
                                     Weight = c(0.1, 0.5, 1, 2, 3, 5, 10)),
      rf = expand.grid(mtry = rf_mtry),
      kknn = expand.grid(kmax = c(5, 7, 9, 11, 13),
                         distance = 2,
                         kernel = "optimal"),
      adaboost = expand.grid(nIter = c(50, 100, 150, 200, 250),
                             method = c("Adaboost.M1", "Real adaboost")),
      LogitBoost = expand.grid(nIter = c(11, 21, 31, 41, 51, 61, 71, 81, 91, 101))
    ))
  }

  list(
    nb = expand.grid(fL = c(0, 1),
                     usekernel = TRUE,
                     adjust = c(0.75, 1, 1.25)),
    svmRadialWeights = expand.grid(sigma = c(0.001, 0.01, 0.05),
                                   C = c(1, 5, 10),
                                   Weight = c(0.5, 1, 2)),
    rf = expand.grid(mtry = rf_mtry),
    kknn = expand.grid(kmax = c(5, 9, 13),
                       distance = 2,
                       kernel = "optimal"),
    adaboost = expand.grid(nIter = c(50, 100, 150),
                           method = c("Adaboost.M1", "Real adaboost")),
    LogitBoost = expand.grid(nIter = c(21, 51, 81))
  )
}

#' Validate and return category class levels with the positive class first
#' @keywords internal
category_class_levels <- function(positive_class = "Y") {
  if (is.null(positive_class) || length(positive_class) != 1L ||
      is.na(positive_class) || !positive_class %in% c("Y", "N")) {
    stop("positive_class must be either 'Y' or 'N'", call. = FALSE)
  }
  c(positive_class, setdiff(c("Y", "N"), positive_class))
}

#' Fit a train-only preprocessing recipe for category ML data
#' @keywords internal
fit_category_preprocess_recipe <- function(x,
                                           common_feature,
                                           label = "Training data",
                                           positive_class = "Y") {
  class_levels <- category_class_levels(positive_class)
  x <- normalize_ml_data_columns(x, label)
  common_feature <- normalize_ml_feature_names(common_feature)
  missing <- setdiff(common_feature, colnames(x))
  if (length(missing) > 0) {
    stop(paste0(label, " is missing required columns: ",
                paste(missing, collapse = ", ")), call. = FALSE)
  }

  x <- x[, common_feature, drop = FALSE]
  assert_unique_ids(x, label)
  x$Var <- factor(as.character(x$Var), levels = class_levels)
  if (any(is.na(x$Var))) {
    stop(paste0(label, " Var must contain only Y/N values"), call. = FALSE)
  }

  feature_names <- common_feature[-c(1:2)]
  x <- coerce_numeric_columns(x, feature_names, label)
  impute_values <- vapply(feature_names, function(feature) {
    col <- x[[feature]]
    if (all(is.na(col))) {
      stop(paste0(label, " has an all-NA feature column: ", feature),
           call. = FALSE)
    }
    mean(col, na.rm = TRUE)
  }, numeric(1))

  na_count <- sum(is.na(x[, feature_names, drop = FALSE]))
  if (na_count > 0) {
    warning(paste0(label, ": ", na_count,
                   " NA values imputed with column means"),
            call. = FALSE)
    for (feature in feature_names) {
      missing_idx <- is.na(x[[feature]])
      if (any(missing_idx)) {
        x[[feature]][missing_idx] <- impute_values[[feature]]
      }
    }
  }
  rownames(x) <- x$ID

  recipe <- structure(
    list(
      common_feature = common_feature,
      feature_names = feature_names,
      impute_values = impute_values,
      positive_class = positive_class,
      class_levels = class_levels
    ),
    class = c("ikl_category_preprocess_recipe", "list")
  )
  list(data = x, recipe = recipe)
}

#' Apply a train-only category ML preprocessing recipe
#' @keywords internal
apply_category_preprocess_recipe <- function(x,
                                             recipe,
                                             label = "data",
                                             common_feature = NULL) {
  if (is.null(common_feature)) {
    common_feature <- recipe$common_feature
  }
  x <- normalize_ml_data_columns(x, label)
  common_feature <- normalize_ml_feature_names(common_feature)
  missing <- setdiff(common_feature, colnames(x))
  if (length(missing) > 0) {
    stop(paste0(label, " is missing required columns: ",
                paste(missing, collapse = ", ")), call. = FALSE)
  }

  x <- x[, common_feature, drop = FALSE]
  assert_unique_ids(x, label)
  x$Var <- factor(as.character(x$Var), levels = recipe$class_levels)
  if (any(is.na(x$Var))) {
    stop(paste0(label, " Var must contain only Y/N values"), call. = FALSE)
  }

  feature_names <- recipe$feature_names
  x <- coerce_numeric_columns(x, feature_names, label)
  na_count <- sum(is.na(x[, feature_names, drop = FALSE]))
  if (na_count > 0) {
    warning(paste0(label, ": ", na_count,
                   " NA values imputed with column means from training preprocessing recipe"),
            call. = FALSE)
    for (feature in feature_names) {
      missing_idx <- is.na(x[[feature]])
      if (any(missing_idx)) {
        x[[feature]][missing_idx] <- recipe$impute_values[[feature]]
      }
    }
  }
  rownames(x) <- x$ID
  x
}

#' IMPRES helper runtime dependencies
#' @keywords internal
impres_required_packages <- function() {
  c(
    "matrixStats", "plyr", "ppcor", "survival", "ROCR", "Hmisc",
    "rms", "mixtools", "lme4", "lmerTest", "plotrix"
  )
}

#' Load package exports into an isolated IMPRES source environment
#' @keywords internal
load_impres_packages_into_env <- function(env, packages = impres_required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(paste0(
      "IMPRES helpers require missing packages: ",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }

  for (pkg in packages) {
    for (nm in getNamespaceExports(pkg)) {
      assign(nm, getExportedValue(pkg, nm), envir = env)
    }
  }
  invisible(env)
}

#' Source packaged IMPRES helper scripts without attaching packages globally
#' @keywords internal
source_impres_extdata_files <- function(env) {
  files <- c(
    "ImmRes_output.R",
    "ImmRes_generic.R",
    "ImmRes_OE.R",
    "ImmRes_master.R",
    "ImmRes1_denovoCellTypeSig.R",
    "ImmRes2_immuneResistanceProgram.R",
    "ImmRes3_longitudinal.R",
    "ImmRes4_predictICBresponses.R",
    "ImmRes5_valCohort2.R",
    "ImmRes6_pancanDrug.R"
  )

  for (file in files) {
    sys.source(system.file("extdata", file, package = "iklSurvML", mustWork = TRUE),
               envir = env)
  }
  invisible(env)
}

#' Load IMPRES helpers into a private environment
#' @keywords internal
load_impres_extdata_helpers <- function(parent = globalenv()) {
  env <- new.env(parent = parent)
  load_impres_packages_into_env(env)
  source_impres_extdata_files(env)
  env
}
