# Debug C-index differences
library(testthat)

source("../R/ml-algorithms.R")
source("../R/helpers.R")
source("../R/ml-unicox.R")
source("../R/ml-combinations.R")
source("../R/ml-dev-prog-sig.R")
source("../R/ml-dev-prog-sig-optimized.R")
source("../../Mime-original/R/ML.Dev.Prog.Sig.R")

set.seed(5201314)
n_samples <- 100
n_genes <- 10

train_data <- data.frame(
  ID = paste0("Sample", 1:n_samples),
  OS.time = rexp(n_samples, rate = 0.01) + 100,
  OS = rbinom(n_samples, 1, 0.6)
)
for (i in 1:n_genes) {
  train_data[[paste0("Gene", i)]] <- rnorm(n_samples, mean = 5, sd = 2) +
    0.3 * (train_data$OS.time - mean(train_data$OS.time)) / sd(train_data$OS.time)
}

set.seed(5201315)
val_data <- data.frame(
  ID = paste0("Val", 1:50),
  OS.time = rexp(50, rate = 0.01) + 100,
  OS = rbinom(50, 1, 0.5)
)
for (i in 1:n_genes) {
  val_data[[paste0("Gene", i)]] <- rnorm(50, mean = 5, sd = 2)
}

list_train_vali_Data <- list(validation = val_data)
candidate_genes <- paste0("Gene", 1:n_genes)

cat("Testing single CoxBoost model...\n")

# Run original - single CoxBoost
cat("\n--- Original CoxBoost ---\n")
result_orig <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = candidate_genes,
  unicox.filter.for.candi = FALSE,
  mode = "single",
  single_ml = "CoxBoost",
  seed = 5201314,
  nodesize = 5,
  cores_for_parallel = 1
)
cat("Original C-index:", result_orig$Cindex.res$Cindex, "\n")
cat("Original RS values (first 5):", head(result_orig$riskscore[[1]]$validation$RS, 5), "\n")

# Run optimized - single CoxBoost
cat("\n--- Optimized CoxBoost ---\n")
result_opt <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = candidate_genes,
  unicox.filter.for.candi = FALSE,
  mode = "single",
  single_ml = "CoxBoost",
  seed = 5201314,
  nodesize = 5,
  cores_for_parallel = 1
)
cat("Optimized C-index:", result_opt$Cindex.res$Cindex, "\n")
cat("Optimized RS values (first 5):", head(result_opt$riskscore[[1]]$validation$RS, 5), "\n")

# Compare
cat("\n--- Comparison ---\n")
rs_orig <- result_orig$riskscore[[1]]$validation$RS
rs_opt <- result_opt$riskscore[[1]]$validation$RS
cat("RS correlation:", cor(rs_orig, rs_opt), "\n")
cat("RS max diff:", max(abs(rs_orig - rs_opt)), "\n")
cat("C-index diff:", abs(result_orig$Cindex.res$Cindex - result_opt$Cindex.res$Cindex), "\n")
