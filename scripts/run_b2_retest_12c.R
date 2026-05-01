artifact_dir <- Sys.getenv("IKLSURVML_B2_ARTIFACT_DIR", "/tmp/iklSurvML-retest/B2-12c-lab")
workspace <- normalizePath(Sys.getenv("IKLSURVML_B2_WORKSPACE", getwd()), mustWork = TRUE)
dataset_path <- Sys.getenv("IKLSURVML_B2_FIXTURE", "")

dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

warnings_seen <- character()
error_message <- NULL
result <- NULL
elapsed_seconds <- NA_real_

writeLines(character(), file.path(artifact_dir, "warnings.txt"))
writeLines(character(), file.path(artifact_dir, "error.txt"))

if (!nzchar(dataset_path)) {
  stop("Set IKLSURVML_B2_FIXTURE to the shared smoke RDS fixture path")
}
dataset <- readRDS(normalizePath(dataset_path, mustWork = TRUE))

pkgload::load_all(workspace, export_all = FALSE, quiet = TRUE)

start_time <- Sys.time()
result <- tryCatch(
  withCallingHandlers(
    ML.Dev.Prog.Sig.Fast(
      train_data = dataset$train_data,
      list_train_vali_Data = dataset$list_train_vali_Data,
      candidate_genes = dataset$candidate_genes,
      unicox.filter.for.candi = FALSE,
      mode = "all",
      seed = dataset$seed,
      nodesize = dataset$nodesize,
      use_parallel = TRUE,
      cores_for_parallel = 12
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
    }
  ),
  error = function(e) {
    error_message <<- conditionMessage(e)
    NULL
  }
)
end_time <- Sys.time()
elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (length(warnings_seen) > 0) {
  writeLines(unique(warnings_seen), file.path(artifact_dir, "warnings.txt"))
}

status <- if (!is.null(error_message)) "errored" else "completed"
writeLines(status, file.path(artifact_dir, "status.txt"))
writeLines(sprintf("%.3f", elapsed_seconds), file.path(artifact_dir, "runtime_seconds.txt"))

if (!is.null(error_message)) {
  writeLines(error_message, file.path(artifact_dir, "error.txt"))
  quit(status = 1)
}

saveRDS(result, file.path(artifact_dir, "result.rds"))

model_names <- names(result$ml.res %||% list())
cindex_rows <- if (!is.null(result$Cindex.res)) nrow(result$Cindex.res) else NA_integer_
object_names <- names(result)
shape_lines <- capture.output(str(result, max.level = 2, list.len = 20))

writeLines(
  c(
    sprintf("model_count: %d", length(model_names)),
    sprintf("cindex_row_count: %s", as.character(cindex_rows)),
    sprintf("first_models: %s", paste(utils::head(model_names, 5), collapse = " | ")),
    sprintf("last_models: %s", paste(utils::tail(model_names, 5), collapse = " | ")),
    sprintf("top_level_names: %s", paste(object_names, collapse = ", "))
  ),
  file.path(artifact_dir, "inspection.txt")
)
writeLines(shape_lines, file.path(artifact_dir, "result-structure.txt"))

cat(sprintf("STATUS=%s\n", status))
cat(sprintf("RUNTIME_SEC=%.3f\n", elapsed_seconds))
cat(sprintf("MODEL_COUNT=%d\n", length(model_names)))
cat(sprintf("CINDEX_ROWS=%s\n", as.character(cindex_rows)))
