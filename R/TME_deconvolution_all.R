#' Tumor microenvironment deconvolution using immunedeconv
#'
#' @param inputmatrix.list A gene expression dataframe. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead.
#' @param deconvolution_method Deconvolution Methods in immunedeconv::deconvolution_methods
#' @param microarray_names Please tell us which datasets are microarray, use the names of elements in inputmatrix.list. such as c("CGGA.array", "GSE108474", "GSE16011", "GSE43289", "GSE7696") if none, enter "none".
#' @param column Only relevant if `gene_expression` is an ExpressionSet. Defines in which column
#'   of fData the HGNC symbol can be found.
#' @param indications a character vector with one indication per
#'   sample for TIMER and ConsensusTME(immunedeconv::timer_available_cancers).
#'   Argument is ignored for all other methods.
#' @param tumor use a signature matrix/procedure optimized for tumor samples,
#'   if supported by the method. Currently affects EPIC and
#' @param rmgenes a character vector of gene symbols. Exclude these genes from the analysis.
#'   Use this to exclude e.g. noisy genes.
#' @param scale_mrna logical. If FALSE, disable correction for mRNA content of different cell types.
#'   This is supported by methods that compute an absolute score (EPIC)
#' @param expected_cell_types Limit the analysis to the cell types given in this list. If the cell
#'   types present in the sample are known *a priori*, setting this can improve results for
#'   xCell (see https://github.com/grst/immunedeconv/issues/1).
#' @param ... arguments passed to the respective method
#' @return `list` containing `data.frame` with `cell_type` as first column and a column with the
#'     calculated cell fractions for each sample.
#' @return a list containing deconvolution scores in each cohorts and ML methods
#' @export
#'
#' @examples
#' \dontrun{
#' test.devo <- TME_deconvolution_all(list_train_vali_Data)
#' }
TME_deconvolution_all <- function(inputmatrix.list, # A list contain the dataframes (colnames:ID,OS.time,OS,other genes), log2(x+1)， OS.time(day), OS(0/1)
                                  deconvolution_method = c("xcell", "epic", "abis", "estimate", "cibersort", "cibersort_abs"), # Deconvolution Methods in c("quantiseq", "xcell", "epic", "abis", "mcp_counter", "estimate", "cibersort", "cibersort_abs", "timer", "consensus_tme")
                                  microarray_names = "none", # Please tell us which datasets are microarray, use the names of elements in inputmatrix.list. such as c("CGGA.array", "GSE108474", "GSE16011", "GSE43289", "GSE7696") if none, enter "none".
                                  indications = NULL,
                                  tumor = TRUE,
                                  column = "gene_symbol",
                                  rmgenes = NULL,
                                  scale_mrna = TRUE,
                                  expected_cell_types = NULL,
                                  ...) {
  message("--- Data preprocessing ---")
  annotate_cell_type <- function(result_table, method) {
    cell_type_map %>%
      filter(method_dataset == !!method) %>%
      inner_join(result_table, by = "method_cell_type") %>%
      dplyr::select(-method_cell_type, -method_dataset)
  }
  immunedeconv::set_cibersort_binary(system.file("extdata", "CIBERSORT.R", package = "iklSurvML"))
  immunedeconv::set_cibersort_mat(system.file("extdata", "LM22.txt", package = "iklSurvML"))
  # Deconvolution by immunedeconv
  TME_deconvolution <- function(gene_expression,
                                deconvolution_method = c("xcell", "epic", "abis", "estimate", "cibersort", "cibersort_abs"),
                                tumor = TRUE,
                                arrays = FALSE, column = "gene_symbol",
                                rmgenes = NULL, scale_mrna = TRUE,
                                expected_cell_types = NULL,
                                ...) {
    deconvolution_methods <- c(
      "quantiseq", "xcell", "epic", "abis", "mcp_counter", "estimate", "cibersort", "cibersort_abs", "timer", "consensus_tme"
    )

    if (all(deconvolution_method %in% deconvolution_methods)) {
      # convert expression set to matrix, if required.
      if (is(gene_expression, "ExpressionSet")) {
        gene_expression <- immunedeconv::eset_to_matrix(gene_expression, column)
      }

      if (!is.null(rmgenes)) {
        gene_expression <- gene_expression[!rownames(gene_expression) %in% rmgenes, ]
      }

      tme_combine <- list()

      for (method in deconvolution_method) {
        message(paste0("\n", ">>> Running ", method))

        # run ed method
        res <- switch(method,
          xcell = immunedeconv::deconvolute_xcell(gene_expression, arrays = arrays, expected_cell_types = expected_cell_types, ...),
          epic = immunedeconv::deconvolute_epic(gene_expression, tumor = tumor, scale_mrna = scale_mrna, ...),
          cibersort = immunedeconv::deconvolute_cibersort(gene_expression,
            absolute = FALSE,
            arrays = arrays, ...
          ),
          cibersort_abs = immunedeconv::deconvolute_cibersort(gene_expression,
            absolute = TRUE,
            arrays = arrays, ...
          ),
          abis = immunedeconv::deconvolute_abis(gene_expression, arrays = arrays),
          estimate = immunedeconv::deconvolute_estimate(gene_expression),
          timer = immunedeconv::deconvolute(gene_expression, "timer",
                                            indications=indications),
          consensus_tme = immunedeconv::deconvolute(gene_expression, "consensus_tme",
                                                    indications=indications),
          quantiseq = immunedeconv::deconvolute(gene_expression, "quantiseq"),
          mcp_counter = immunedeconv::deconvolute(gene_expression, "mcp_counter")
        )

        # convert to tibble and annotate unified cell_type names
        res <- res %>%
          as_tibble(rownames = "method_cell_type") %>%
          annotate_cell_type(method = method)

        tme_combine[[method]] <- res
      }

      resultList <- list("tme_combine" = tme_combine)
      return(resultList)
    } else {
      invalid_methods <- setdiff(deconvolution_method, deconvolution_methods)
      stop(paste0(
        "Unsupported deconvolution_method value(s): ",
        paste(invalid_methods, collapse = ", ")
      ), call. = FALSE)
    }
  }

  inputmatrix.list <- lapply(inputmatrix.list, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)
    })
    return(x)
  })
  tme_decon_list <- list()
  # cat("Datasets names:", paste(names(inputmatrix.list), collapse = ", "), "\n")
  # cat("Please tell us which datasets are microarray, if none, skip the step\n")  # 提供一个示例


  selected_columns <- microarray_names
  if (!is.character(selected_columns) || length(selected_columns) == 0L ||
      any(is.na(selected_columns))) {
    stop("microarray_names must be 'none' or one or more dataset names", call. = FALSE)
  }
  microarray_is_none <- length(selected_columns) == 1L && identical(selected_columns, "none")
  if (!microarray_is_none) {
    invalid_microarray_names <- setdiff(selected_columns, names(inputmatrix.list))
    if (length(invalid_microarray_names) > 0L) {
      stop(paste0(
        "microarray_names contains dataset names not present in inputmatrix.list: ",
        paste(invalid_microarray_names, collapse = ", ")
      ), call. = FALSE)
    }
  }

  tme_errors <- character()
  for (i in seq_along(inputmatrix.list)) {
    dataset_name <- names(inputmatrix.list)[i]
    train_data <- inputmatrix.list[[i]]
    test.matrix <- train_data %>%
      magrittr::set_rownames(.$ID) %>%
      dplyr::select(-c(ID, OS.time, OS)) %>%
      t()

    resultList <- tryCatch(
      {
        use_microarray <- !microarray_is_none && dataset_name %in% selected_columns
        TME_deconvolution(
          gene_expression = test.matrix,
          deconvolution_method = deconvolution_method,
          arrays = use_microarray,
          tumor = tumor,
          column = column,
          rmgenes = rmgenes,
          scale_mrna = scale_mrna,
          expected_cell_types = expected_cell_types,
          ...
        )
      },
      error = function(e) {
        tme_errors <<- c(tme_errors, paste0(dataset_name, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(resultList)) {
      print("Success")
      tme_decon_list[[dataset_name]] <- resultList
    }
  }

  if (length(tme_errors) > 0L) {
    attr(tme_decon_list, "errors") <- tme_errors
    stop(paste0(
      "TME deconvolution failed for ",
      length(tme_errors),
      " cohort(s): ",
      paste(tme_errors, collapse = "; ")
    ), call. = FALSE)
  }
  return(tme_decon_list)
}
