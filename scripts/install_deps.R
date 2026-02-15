# Install dependencies for Mime package
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

# List of required packages
required_packages <- c(
  "CoxBoost",
  "randomForestSRC",
  "glmnet",
  "plsRcox",
  "superpc",
  "gbm",
  "survivalsvm",
  "survival",
  "testthat"
)

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install each package
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch({
      if (pkg %in% c("CoxBoost", "randomForestSRC", "plsRcox")) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg)
      }
    }, error = function(e) {
      cat(sprintf("Failed to install %s: %s\n", pkg, e$message))
    })
  } else {
    cat(sprintf("%s already installed\n", pkg))
  }
}

cat("\nAll dependencies installed!\n")
