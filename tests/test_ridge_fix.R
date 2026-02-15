# Quick test for Ridge fix
source("../R/ml-algorithms.R")
source("../R/helpers.R")
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

candidate_genes <- paste0("Gene", 1:n_genes)
pre_var <- candidate_genes
est_dd <- train_data[, c("OS.time", "OS", candidate_genes)]
val_dd_list <- list(validation = val_data[, c("OS.time", "OS", candidate_genes)])

cat("Testing Ridge directly...\n\n")

# Original Ridge implementation
cat("--- Original Ridge ---\n")
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(survival::Surv(est_dd$OS.time, est_dd$OS))
set.seed(5201314)
fit_orig <- glmnet::glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
cv_orig <- glmnet::cv.glmnet(x1, x2, nfold = 10, family = "cox")

rs_orig <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit_orig, type = "response",
                                          newx = as.matrix(x[, -c(1,2)]),
                                          s = cv_orig$lambda.min)))
})
cindex_orig <- as.numeric(summary(survival::coxph(survival::Surv(OS.time, OS) ~ RS, rs_orig$validation))$concordance[1])
cat("Original C-index:", cindex_orig, "\n")
cat("Original lambda.min:", cv_orig$lambda.min, "\n")

# Optimized Ridge implementation
cat("\n--- Optimized Ridge ---\n")
fit_opt <- train_ridge(est_dd, pre_var, seed = 5201314)

rs_opt <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = predict_ridge(fit_opt, x, pre_var))
})
cindex_opt <- as.numeric(summary(survival::coxph(survival::Surv(OS.time, OS) ~ RS, rs_opt$validation))$concordance[1])
cat("Optimized C-index:", cindex_opt, "\n")
cat("Optimized lambda.min:", fit_opt$cv.fit$lambda.min, "\n")

cat("\n--- Comparison ---\n")
cat("C-index difference:", abs(cindex_orig - cindex_opt), "\n")
cat("lambda.min match:", cv_orig$lambda.min == fit_opt$cv.fit$lambda.min, "\n")
cat("RS correlation:", cor(rs_orig$validation$RS, rs_opt$validation$RS), "\n")

if (abs(cindex_orig - cindex_opt) < 0.0001) {
  cat("\n✓ Ridge fix successful!\n")
} else {
  cat("\n✗ Ridge still has differences\n")
}
