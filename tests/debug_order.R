# Debug: Compare specific models in all mode
source("../R/ml-algorithms.R")
source("../R/helpers.R")
source("../R/ml-unicox.R")
source("../R/ml-combinations.R")
source("../R/ml-dev-prog-sig.R")
source("../R/ml-dev-prog-sig-optimized.R")
source("../../Mime-original/R/ML.Dev.Prog.Sig.R")

set.seed(5201314)
n_samples <- 80
n_genes <- 8

train_data <- data.frame(
  ID = paste0("Sample", 1:n_samples),
  OS.time = rexp(n_samples, rate = 0.01) + 100,
  OS = rbinom(n_samples, 1, 0.6)
)
for (i in 1:n_genes) {
  train_data[[paste0("Gene", i)]] <- rnorm(n_samples, mean = 5, sd = 2)
}

set.seed(5201315)
val_data <- data.frame(
  ID = paste0("Val", 1:40),
  OS.time = rexp(40, rate = 0.01) + 100,
  OS = rbinom(40, 1, 0.5)
)
for (i in 1:n_genes) {
  val_data[[paste0("Gene", i)]] <- rnorm(40, mean = 5, sd = 2)
}

list_train_vali_Data <- list(validation = val_data)
candidate_genes <- paste0("Gene", 1:n_genes)

cat("Running Original (all mode)...\n")
r_orig <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = candidate_genes,
  unicox.filter.for.candi = FALSE,
  mode = "all",
  seed = 5201314,
  nodesize = 5,
  cores_for_parallel = 1
)

cat("\nRunning Optimized (all mode)...\n")
r_opt <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = candidate_genes,
  unicox.filter.for.candi = FALSE,
  mode = "all",
  seed = 5201314,
  nodesize = 5,
  cores_for_parallel = 1
)

cat("\n========================================\n")
cat("Comparing all models\n")
cat("========================================\n")

orig_names <- names(r_orig$ml.res)
opt_names <- names(r_opt$ml.res)
common <- intersect(orig_names, opt_names)

cat("Original models:", length(orig_names), "\n")
cat("Optimized models:", length(opt_names), "\n")
cat("Common models:", length(common), "\n")

# Compare C-index for each model
cat("\n--- C-index differences per model ---\n")
diffs <- data.frame()
for (model in common) {
  orig_ci <- r_orig$Cindex.res$Cindex[r_orig$Cindex.res$Model == model]
  opt_ci <- r_opt$Cindex.res$Cindex[r_opt$Cindex.res$Model == model]

  if (length(orig_ci) > 0 && length(opt_ci) > 0) {
    diff <- abs(orig_ci[1] - opt_ci[1])
    if (diff > 0.001) {
      cat(sprintf("%s: orig=%.6f, opt=%.6f, diff=%.6f\n", model, orig_ci[1], opt_ci[1], diff))
    }
    diffs <- rbind(diffs, data.frame(model = model, diff = diff))
  }
}

cat("\n--- Summary ---\n")
cat("Models with diff > 0.001:", sum(diffs$diff > 0.001), "/", nrow(diffs), "\n")
cat("Max diff:", max(diffs$diff), "\n")
cat("Mean diff:", mean(diffs$diff), "\n")
