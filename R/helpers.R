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
preprocess_data_list <- function(list_train_vali_Data, common_feature) {
  # Select common features
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, common_feature]
  })

  # Convert gene columns to numeric
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, as.numeric)
    return(x)
  })


  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x[, c(2:3)] <- apply(x[, c(2:3)], 2, as.numeric)
    return(x)
  })

  # Remove NA in OS.time and OS
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x <- x[!is.na(x$OS.time) & !is.na(x$OS), ]
    return(x)
  })

  # Remove OS.time <= 0
  list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
    x <- x[x$OS.time > 0, ]
    return(x)
  })

  # Replace NA with mean (with warning)
  for (nm in names(list_train_vali_Data)) {
    x <- list_train_vali_Data[[nm]]
    na_count <- sum(is.na(x[, -c(1:3)]))
    if (na_count > 0) {
      warning(paste0("Dataset '", nm, "': ", na_count,
                     " NA values imputed with column means"))
    }
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(col) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      return(col)
    })
    list_train_vali_Data[[nm]] <- x
  }

  return(list_train_vali_Data)
}

#' Preprocess training data for ML analysis
#'
#' @param train_data Training data
#' @param common_feature Common features to select
#' @return Preprocessed training data
#' @keywords internal
preprocess_train_data <- function(train_data, common_feature) {
  train_data <- train_data[, common_feature]
  train_data[, -c(1:3)] <- apply(train_data[, -c(1:3)], 2, as.numeric)
  train_data[, c(2:3)] <- apply(train_data[, c(2:3)], 2, as.numeric)
  return(train_data)
}
