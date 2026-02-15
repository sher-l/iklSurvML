# Test consistency between original and optimized versions
# Compares actual results, not just counts

library(testthat)

# Source refactored files
source("../R/ml-algorithms.R")
source("../R/helpers.R")
source("../R/ml-unicox.R")
source("../R/ml-combinations.R")
source("../R/ml-dev-prog-sig.R")
source("../R/ml-dev-prog-sig-optimized.R")

# Source original files
source("../../Mime-original/R/ML.Dev.Prog.Sig.R")

# Create test data with known seed for reproducibility
set.seed(5201314)
n_samples <- 150
n_genes <- 15

# Training data
train_data <- data.frame(
  ID = paste0("Sample", 1:n_samples),
  OS.time = rexp(n_samples, rate = 0.01) + 100,
  OS = rbinom(n_samples, 1, 0.6)
)

# Create features with some correlation to survival
for (i in 1:n_genes) {
  if (i <= 5) {
    train_data[[paste0("Gene", i)]] <- rnorm(n_samples, mean = 5, sd = 2) +
      0.3 * (train_data$OS.time - mean(train_data$OS.time)) / sd(train_data$OS.time)
  } else {
    train_data[[paste0("Gene", i)]] <- rnorm(n_samples, mean = 5, sd = 2)
  }
}

# Validation data
set.seed(5201315)
val_data <- data.frame(
  ID = paste0("Val", 1:80),
  OS.time = rexp(80, rate = 0.01) + 100,
  OS = rbinom(80, 1, 0.5)
)
for (i in 1:n_genes) {
  val_data[[paste0("Gene", i)]] <- rnorm(80, mean = 5, sd = 2)
}

list_train_vali_Data <- list(validation = val_data)
candidate_genes <- paste0("Gene", 1:n_genes)

cat("\n========================================\n")
cat("Comparing Original vs Optimized Results\n")
cat("========================================\n\n")

# Run original version
cat("Running Original version...\n")
start_orig <- Sys.time()
result_orig <- tryCatch({
  ML.Dev.Prog.Sig(
    train_data = train_data,
    list_train_vali_Data = list_train_vali_Data,
    candidate_genes = candidate_genes,
    unicox.filter.for.candi = FALSE,
    mode = "all",
    seed = 5201314,
    nodesize = 5,
    cores_for_parallel = 1
  )
}, error = function(e) {
  cat("Error in Original version:", conditionMessage(e), "\n")
  return(NULL)
})
time_orig <- as.numeric(difftime(Sys.time(), start_orig, units = "secs"))
cat(sprintf("Original completed in %.2f seconds\n", time_orig))

# Run optimized version
cat("\nRunning Optimized version...\n")
start_opt <- Sys.time()
result_opt <- tryCatch({
  ML.Dev.Prog.Sig.Fast(
    train_data = train_data,
    list_train_vali_Data = list_train_vali_Data,
    candidate_genes = candidate_genes,
    unicox.filter.for.candi = FALSE,
    mode = "all",
    seed = 5201314,
    nodesize = 5,
    cores_for_parallel = 1
  )
}, error = function(e) {
  cat("Error in Optimized version:", conditionMessage(e), "\n")
  return(NULL)
})
time_opt <- as.numeric(difftime(Sys.time(), start_opt, units = "secs"))
cat(sprintf("Optimized completed in %.2f seconds\n", time_opt))

if (is.null(result_orig) || is.null(result_opt)) {
  cat("\nOne or both versions failed\n")
  quit(save = "no", status = 1)
}

cat("\n========================================\n")
cat("Comparison Results\n")
cat("========================================\n\n")

# Compare model counts
cat("--- Model Counts ---\n")
cat(sprintf("Original models: %d\n", length(result_orig$ml.res)))
cat(sprintf("Optimized models: %d\n", length(result_opt$ml.res)))

# Get common model names
orig_names <- names(result_orig$ml.res)
opt_names <- names(result_opt$ml.res)
common_names <- intersect(orig_names, opt_names)
only_orig <- setdiff(orig_names, opt_names)
only_opt <- setdiff(opt_names, orig_names)

cat(sprintf("\nCommon models: %d\n", length(common_names)))
if (length(only_orig) > 0) {
  cat(sprintf("Only in Original: %d (%s)\n", length(only_orig), paste(head(only_orig, 5), collapse = ", ")))
}
if (length(only_opt) > 0) {
  cat(sprintf("Only in Optimized: %d (%s)\n", length(only_opt), paste(head(only_opt, 5), collapse = ", ")))
}

# Compare C-index results
cat("\n--- C-index Comparison ---\n")
if ("Cindex.res" %in% names(result_orig) && "Cindex.res" %in% names(result_opt)) {
  orig_cindex <- result_orig$Cindex.res
  opt_cindex <- result_opt$Cindex.res

  # Compare C-index values for common models
  cindex_diff <- c()
  for (model in common_names) {
    orig_val <- orig_cindex$Cindex[orig_cindex$Model == model]
    opt_val <- opt_cindex$Cindex[opt_cindex$Model == model]

    if (length(orig_val) > 0 && length(opt_val) > 0) {
      diff <- abs(orig_val[1] - opt_val[1])
      cindex_diff <- c(cindex_diff, diff)
    }
  }

  if (length(cindex_diff) > 0) {
    cat(sprintf("Models compared: %d\n", length(cindex_diff)))
    cat(sprintf("Max C-index difference: %.6f\n", max(cindex_diff)))
    cat(sprintf("Mean C-index difference: %.6f\n", mean(cindex_diff)))

    # Count matches (difference < 0.001)
    matches <- sum(cindex_diff < 0.001)
    cat(sprintf("C-index matches (diff < 0.001): %d / %d (%.1f%%)\n",
                matches, length(cindex_diff), 100 * matches / length(cindex_diff)))
  }
}

# Compare risk scores for a few models
cat("\n--- Risk Score Comparison (sample models) ---\n")
sample_models <- head(common_names, 5)
for (model in sample_models) {
  if (model %in% names(result_orig$riskscore) && model %in% names(result_opt$riskscore)) {
    orig_rs <- result_orig$riskscore[[model]]
    opt_rs <- result_opt$riskscore[[model]]

    if ("validation" %in% names(orig_rs) && "validation" %in% names(opt_rs)) {
      orig_vals <- orig_rs$validation$RS
      opt_vals <- opt_rs$validation$RS

      if (length(orig_vals) == length(opt_vals)) {
        cor_val <- cor(orig_vals, opt_vals)
        max_diff <- max(abs(orig_vals - opt_vals))
        cat(sprintf("%s: correlation=%.6f, max_diff=%.6f\n", model, cor_val, max_diff))
      }
    }
  }
}

# Performance comparison
cat("\n--- Performance ---\n")
cat(sprintf("Original: %.2f seconds\n", time_orig))
cat(sprintf("Optimized: %.2f seconds\n", time_opt))
cat(sprintf("Speedup: %.2fx\n", time_orig / time_opt))

cat("\n========================================\n")
cat("Test completed\n")
cat("========================================\n")
