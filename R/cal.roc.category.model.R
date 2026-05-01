#' Calculating the receiver operating characteristic after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig    Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal  A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#'
#' @return  A list containing the receiver operating characteristic of each predictive model.
#' @export
#'
cal.roc.category.model <- function(res.by.ML.Dev.Pred.Category.Sig, ### 函数计算结果
	                                   cohort.for.cal # 队列要求第一列为ID,第二列为分类变量Var, 值为Y或者N, 从第三列开始为基因，表达矩阵经过了log2(x+1)处理
) {
  sig <- res.by.ML.Dev.Pred.Category.Sig$sig.gene
  positive_class <- res.by.ML.Dev.Pred.Category.Sig$positive_class
  if (is.null(positive_class)) {
    positive_class <- "Y"
  }
  class_levels <- category_class_levels(positive_class)
  negative_class <- setdiff(c("Y", "N"), positive_class)

  rownames(cohort.for.cal) <- cohort.for.cal$ID
  cohort.for.cal <- normalize_ml_data_columns(cohort.for.cal, "cohort.for.cal")
  sig <- normalize_ml_feature_names(sig)

  required_columns <- c("ID", "Var", sig)
  missing_columns <- setdiff(required_columns, colnames(cohort.for.cal))
  if (length(missing_columns) > 0L) {
    stop(paste0(
      "cohort.for.cal is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    ), call. = FALSE)
  }
  if (!identical(colnames(cohort.for.cal)[1:2], c("ID", "Var"))) {
    stop("cohort.for.cal first 2 columns must be ID, Var", call. = FALSE)
  }

  if (all(is.element(sig, colnames(cohort.for.cal))) & identical(colnames(cohort.for.cal)[1:2], c("ID", "Var"))) {
    recipe <- res.by.ML.Dev.Pred.Category.Sig$Preprocess.recipe
    if (!is.null(recipe)) {
      validation <- apply_category_preprocess_recipe(
        cohort.for.cal,
        recipe = recipe,
        label = "cohort.for.cal",
        common_feature = c("ID", "Var", sig)
      )
      validation <- validation[, colnames(validation) %in% c("Var", sig)]
    } else {
      validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
      validation$Var <- factor(validation$Var, levels = class_levels)
    }

    ls_model <- res.by.ML.Dev.Pred.Category.Sig$model
    models <- names(ls_model)
    roc <- lapply(1:length(models), function(i) {
      if (!models[i] == "cancerclass") {
        prob <- predict(ls_model[[models[i]]], validation[, -1], type = "prob") #
        pre <- predict(ls_model[[models[i]]], validation[, -1]) #
        test_set <- data.frame(obs = validation$Var, pred = pre, check.names = FALSE)
        for (cls in class_levels) {
          test_set[[cls]] <- prob[, cls]
        }
        roc <- ROCit::rocit(
          score = test_set[[positive_class]],
          class = test_set$obs,
          negref = negative_class
        )
      } else {
        pData <- data.frame(class = validation$Var, sample = rownames(validation), row.names = rownames(validation))
        phenoData <- Biobase::AnnotatedDataFrame(data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- Biobase::ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)

        prediction <- predict(ls_model[[models[i]]], Sig.Exp.test, positive_class, ngenes = nrow(Sig.Exp), dist = "cor")
        roc <- category_cancerclass_roc(
          observed = validation$Var,
          prediction = prediction,
          positive_class = positive_class
        )
      }
    })

    names(roc) <- models


    return(roc)
  } else {
    stop("Please provide the correct parameters and the cohorts with correct format", call. = FALSE)
  }
}
