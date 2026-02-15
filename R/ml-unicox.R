# Univariate Cox Regression for Feature Selection
#
# This file contains functions for univariate Cox regression analysis

#' Perform univariate Cox regression for gene selection
#'
#' @param gene_list A character vector of gene names to analyze
#' @param inputSet A data frame with ID, OS.time, OS, and gene expression columns
#' @param unicox_pcutoff P-value threshold for selecting significant genes (default 0.05)
#' @return A character vector of selected gene names
#' @keywords internal
sig_unicox <- function(gene_list, inputSet, unicox_pcutoff = 0.05) {
  message("Starting the data preprocess")

  # Data preprocessing
  message("Rejecting a null value")

  # Replace NA with 0
  inputSet[is.na(inputSet)] <- 0
  inputSet <- as.data.frame(inputSet)
  inputSet$OS.time <- as.numeric(inputSet$OS.time)

  # Remove rows with OS.time <= 0
  inputSet <- inputSet[inputSet$OS.time > 0, ]

  # Standardize gene names
  gene_list <- gsub("-", ".", gene_list)
  gene_list <- gsub("_", ".", gene_list)
  colnames(inputSet)[4:ncol(inputSet)] <- gsub("-",
    ".",
    colnames(inputSet)[4:ncol(inputSet)]
  )
  colnames(inputSet)[4:ncol(inputSet)] <- gsub("_",
    ".",
    colnames(inputSet)[4:ncol(inputSet)]
  )

  message("Gets the intersection of genelist and expression profile")
  # Get intersection of gene list and expression profile
  comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)], gene_list)

  message("Processing the input representation matrix")
  # Process input expression matrix
  inputSet <- inputSet[, c("ID", "OS.time", "OS", comsa1)]
  inputSet[, c(1:2)] <- apply(inputSet[, c(1:2)], 2, as.factor)
  inputSet[, c(2:ncol(inputSet))] <- apply(inputSet[, c(2:ncol(inputSet))], 2, as.numeric)
  inputSet <- as.data.frame(inputSet)

  message("Data preprocessing completed")

  # Univariate Cox regression
  message("Starting the univariable cox regression")

  unicox <- data.frame()
  for (i in seq_len(ncol(inputSet[, 4:ncol(inputSet)]))) {
    display_progress(
      index = i,
      totalN = ncol(inputSet[, 4:ncol(inputSet)])
    )
    gene <- colnames(inputSet[, 4:ncol(inputSet)])[i]
    tmp <- data.frame(
      expr = as.numeric(inputSet[, 4:ncol(inputSet)][, i]),
      futime = inputSet$OS.time,
      fustat = inputSet$OS,
      stringsAsFactors = FALSE
    )
    cox <- survival::coxph(survival::Surv(futime, fustat) ~ expr, data = tmp)
    coxSummary <- summary(cox)
    unicox <- rbind.data.frame(
      unicox,
      data.frame(
        gene = gene,
        HR = as.numeric(coxSummary$coefficients[, "exp(coef)"])[1],
        z = as.numeric(coxSummary$coefficients[, "z"])[1],
        pvalue = as.numeric(coxSummary$coefficients[, "Pr(>|z|)"])[1],
        lower = as.numeric(coxSummary$conf.int[, 3][1]),
        upper = as.numeric(coxSummary$conf.int[, 4][1]),
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }

  message("Finished the univariable cox regression")

  # Select significant genes
  selgene <- unicox[which(unicox$pvalue < unicox_pcutoff), "gene"]
  return(selgene)
}
