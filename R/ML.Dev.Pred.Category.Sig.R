#' Developing the optimal predictive model for the dichotomous variables with machine learning algorithms
#' 
#' A function can be used to develop the predictive model for dichotomous variables with seven machine learning algorithms.
#' 
#' @param train_data The training data with the 'ID' and 'Var' as the first two columns. Starting in the third column are the variables used to construct the model. 'Var' is the target predictor variable for constructing the model. 'Var' contains only Y or N.
#' @param list_train_vali_Data A list containing the training data and the other validation data. All the validation data have the same data form as the training data.
#' @param candidate_genes The candidate variables used for constructing the predictive model.
#' @param methods There are seven algorithms for developing the predictive model including 'nb', 'svmRadialWeights', 'rf', 'kknn', 'adaboost', 'LogitBoost', 'cancerclass'. 'nb':Naive Bayes algorithm. 'svmRadialWeights': Support Vector Machine (SVM). 'rf': Random Forest. 'kknn': K-nearest Neighbors.'adaboost': AdaBoost Classification Trees. 'LogitBoost':Boosted Logistic Regressions. 'cancerclass': Cancerclass. 
#' @param seed The seed you can set as any positive number, for example, 5201314.
#' @param cores_for_parallel The cores you can choose for parallel operation. The default is 12.The bigger the better if the configuration allows it.
#' @param positive_class The class treated as the positive/event class for ROC/AUC. Defaults to "Y".
#' @param feature_alignment How candidate genes are aligned across cohorts.
#'   "strict" (default) requires all candidate genes to be present in every
#'   training/validation cohort. "intersection" preserves the legacy behavior
#'   of training only on genes shared by all cohorts, with a warning for drops.
#'
#' @return A list containing the predictive model, the AUC, the ROC, and the candidate variables, all of which are developed by each single algorithm.
#'
#' @examples
#' Resolve category ML methods and optional engines
#' @keywords internal
resolve_category_methods <- function(methods = NULL) {
  method_packages <- category_method_packages()
  valid_methods <- names(method_packages)

  if (is.null(methods)) {
    methods <- valid_methods[vapply(method_packages, function(pkgs) {
      all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
    }, logical(1))]
    methods <- intersect(methods, valid_methods)
    if (length(methods) == 0L) {
      stop("No category ML engines are available; install at least one supported engine package.", call. = FALSE)
    }
    return(methods)
  }

  invalid <- setdiff(methods, valid_methods)
  if (length(invalid) > 0) {
    stop(paste0("methods contains unsupported values: ", paste(invalid, collapse = ", ")), call. = FALSE)
  }

  missing_by_method <- vapply(methods, function(method) {
    pkgs <- method_packages[[method]]
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
      paste(missing, collapse = ", ")
    } else {
      NA_character_
    }
  }, character(1))
  missing_by_method <- missing_by_method[!is.na(missing_by_method)]
  if (length(missing_by_method) > 0L) {
    stop(paste0(
      "Unavailable category ML engine package(s): ",
      paste(paste0(names(missing_by_method), " requires ", missing_by_method), collapse = "; ")
    ), call. = FALSE)
  }

  methods
}

#' Validate category ML data before creating workers
#' @keywords internal
validate_category_ml_inputs <- function(train_data,
                                        list_train_vali_Data,
                                        candidate_genes,
                                        common_feature) {
  if (!is.data.frame(train_data)) {
    stop("train_data must be a data.frame", call. = FALSE)
  }
  if (!is.list(list_train_vali_Data) || length(list_train_vali_Data) == 0) {
    stop("list_train_vali_Data must be a non-empty list of data.frames", call. = FALSE)
  }
  if (!identical(colnames(train_data)[1:2], c("ID", "Var"))) {
    stop("first 2 columns of train_data must be ID, Var", call. = FALSE)
  }
  bad_datasets <- names(list_train_vali_Data)[
    !vapply(list_train_vali_Data, function(x) {
      is.data.frame(x) && identical(colnames(x)[1:2], c("ID", "Var"))
    }, logical(1))
  ]
  if (length(bad_datasets) > 0) {
    stop(paste0(
      "first 2 columns of each validation dataset must be ID, Var. Invalid: ",
      paste(bad_datasets, collapse = ", ")
    ), call. = FALSE)
  }
  if (length(candidate_genes) < 2) {
    stop("candidate_genes must contain at least 2 features", call. = FALSE)
  }
  if (!identical(common_feature[1:2], c("ID", "Var")) || length(common_feature) <= 3) {
    stop("candidate_genes must overlap train_data and all validation datasets", call. = FALSE)
  }

  all_var_values <- unique(c(
    as.character(train_data$Var),
    unlist(lapply(list_train_vali_Data, function(x) as.character(x$Var)), use.names = FALSE)
  ))
  invalid_var_values <- setdiff(all_var_values[!is.na(all_var_values)], c("Y", "N"))
  if (any(is.na(all_var_values)) || length(invalid_var_values) > 0) {
    stop("Var must contain only non-missing Y/N values", call. = FALSE)
  }

  invisible(TRUE)
}

#' @export
ML.Dev.Pred.Category.Sig <- function(train_data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N，
                                     list_train_vali_Data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N
                                     candidate_genes = NULL,
                                     methods = NULL, # c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
                                     seed = 5201314, # 5201314
                                     cores_for_parallel = 12, #
                                     positive_class = "Y",
                                     feature_alignment = c("strict", "intersection")
) {
  ###### loading the function #######



  class_levels <- category_class_levels(positive_class)
  negative_class <- setdiff(c("Y", "N"), positive_class)
  feature_alignment <- match.arg(feature_alignment)

  model.Dev <- function(training, method, sig) {
    training <- training[, colnames(training) %in% c("Var", sig)]
    training$Var <- factor(as.character(training$Var), levels = class_levels)
    # 7 models adpoted in this study as followings:
    #' nb': navie bayes
    #' svmRadialWeights': Support Vector Machines with Class Weights
    #' rf': random forest
    #' kknn': k-Nearest Neighbors
    #' adaboost':AdaBoost Classification Trees
    #' LogitBoost':Boosted Logistic Regressions
    #' cancerclass': cancerclass

    # Grid search for parameter tuning
    Grid <- category_tune_grid(n_features = length(sig))
    TuneLength <- list(
      nb = nrow(Grid[["nb"]]),
      svmRadialWeights = nrow(Grid[["svmRadialWeights"]]),
      rf = nrow(Grid[["rf"]]),
      kknn = nrow(Grid[["kknn"]]),
      adaboost = nrow(Grid[["adaboost"]]),
      LogitBoost = nrow(Grid[["LogitBoost"]])
    )
    ls_model <- lapply(method, function(m) {
      if (m == "cancerclass") { # cancerclass is not avaliable in caret
        pData <- data.frame(class = training$Var, sample = rownames(training), row.names = rownames(training))
        phenoData <- Biobase::AnnotatedDataFrame(data = pData)
        Sig.Exp <- t(training[, -1])
        Sig.Exp.train <- Biobase::ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
        predictor <- cancerclass::fit(Sig.Exp.train, method = "welch.test")
        model.tune <- predictor
      } else {
        f <- 5 # f folds resampling
        r <- 10 # r repeats
        n <- f * r

        # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
        seeds <- vector(mode = "list", length = n + 1)
        # the number of tuning parameter
        for (i in 1:n) seeds[[i]] <- sample.int(n = 1000, TuneLength[[m]])

        # for the last model
        seeds[[n + 1]] <- sample.int(1000, 1)


        ctrl <- caret::trainControl(
          method = "repeatedcv",
          number = f, ## 5-folds cv
          summaryFunction = caret::twoClassSummary, # Use AUC to pick the best model
          classProbs = TRUE,
          repeats = r, ## 10-repeats cv,
          seeds = seeds
        )



        model.tune <- caret::train(Var ~ .,
          data = training,
          method = m,
          metric = "ROC",
          trControl = ctrl,
          tuneGrid = Grid[[m]]
        )
      }
      print(m)
      return(model.tune)
    })


    names(ls_model) <- method

    return(ls_model)
  }





  # CompareModel <- function(training, validation, method,sig){
  #
  #   training <- training[,colnames(training) %in% c('Var', sig)]
  #   validation  <- validation[,colnames(validation) %in% c('Var',sig)]
  #
  #
  #
  #
  #
  #
  #
  #
  #   auc <- lapply(ls_model,function(model.tune){
  #     if(class(model.tune) == 'predictor'){
  #       pData <- data.frame(class = validation$Var, sample = rownames(validation),row.names = rownames(validation))
  #       phenoData <- new("AnnotatedDataFrame",data=pData)
  #       Sig.Exp <- t(validation[,-1])
  #       Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
  #       prediction <- predict(model.tune, Sig.Exp.test,"N", ngenes=nrow(Sig.Exp), dist = "cor")
  #       roc <- roc(response  = prediction@prediction[,'class_membership'],
  #                  predictor = as.numeric(prediction@prediction[,'z'])
  #       )
  #       roc_result <- coords(roc, "best")
  #       auc <- data.frame(ROC=as.numeric(roc$auc), Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
  #     }else {
  #       prob <- predict(model.tune,validation[,-1],type = "prob")
  #       pre <- predict(model.tune,validation[,-1])
  #       test_set <- data.frame(obs = validation$Var, N = prob[,'N'], Y = prob[,'Y'], pred=pre)
  #       auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
  #     }
  #
  #     return(auc)
  #   }) %>% base::do.call(rbind,.)
  #
  #   rownames(auc) <- method
  #
  #   res <- list()
  #
  #   names(ls_model) = method
  #
  #   res[['model']] <- ls_model
  #   res[['auc']] <- auc
  #
  #
  #   return(res)
  #
  # }


  cal.model.auc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = class_levels)
    
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    auc <- lapply(1:length(models), function(i) {
      if (models[i] == "cancerclass") {
        model.tune <- ls_model[[i]]
        pData <- data.frame(class = validation$Var, sample = rownames(validation), row.names = rownames(validation))
        phenoData <- Biobase::AnnotatedDataFrame(data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- Biobase::ExpressionSet(assayData = as.matrix(Sig.Exp), phenoData = phenoData)
        prediction <- predict(model.tune, Sig.Exp.test, positive_class, ngenes = nrow(Sig.Exp), dist = "cor")
        roc <- pROC::roc(
          response = prediction@prediction[, "class_membership"],
          predictor = as.numeric(prediction@prediction[, "z"])
        )
        roc_result <- pROC::coords(roc, "best")
        auc <- data.frame(ROC = as.numeric(roc$auc), Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
      } else {
        model.tune <- ls_model[[i]]
        prob <- predict(model.tune, validation[, -1], type = "prob")
        pre <- predict(model.tune, validation[, -1])
        test_set <- data.frame(obs = validation$Var, pred = pre, check.names = FALSE)
        for (cls in class_levels) {
          test_set[[cls]] <- prob[, cls]
        }
        auc <- caret::twoClassSummary(test_set, lev = class_levels)
      }
      
      return(auc)
    }) %>% base::do.call(rbind, .)
    
    rownames(auc) <- names(ls_model)
    
    return(auc)
  }


  cal.model.roc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = class_levels)

    ls_model <- res.by.model.Dev
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
        roc <- pROC::roc(
          response = prediction@prediction[, "class_membership"],
          predictor = as.numeric(prediction@prediction[, "z"])
        )
      }
    })

    names(roc) <- models


    return(roc)
  }



  message("---loading the function---")


  
  data_names <- names(list_train_vali_Data)
  if (is.null(data_names)) {
    data_names <- as.character(seq_along(list_train_vali_Data))
  }
  list_train_vali_Data <- lapply(data_names, function(nm) {
    normalize_ml_data_columns(list_train_vali_Data[[nm]], paste0("Dataset '", nm, "'"))
  })
  names(list_train_vali_Data) <- data_names
  train_data <- normalize_ml_data_columns(train_data, "Training data")
  candidate_genes <- normalize_ml_feature_names(candidate_genes)
  common_feature <- c("ID", "Var", candidate_genes)
  if (identical(feature_alignment, "strict")) {
    datasets <- c(list(`Training data` = train_data), list_train_vali_Data)
    missing_by_dataset <- vapply(names(datasets), function(nm) {
      paste(setdiff(common_feature, colnames(datasets[[nm]])), collapse = ", ")
    }, character(1))
    missing_by_dataset <- missing_by_dataset[nzchar(missing_by_dataset)]
    if (length(missing_by_dataset) > 0L) {
      stop(paste0(
        "feature_alignment='strict' requires ID, Var, and all candidate_genes ",
        "in every training/validation dataset. Missing: ",
        paste(paste0(names(missing_by_dataset), " [", missing_by_dataset, "]"), collapse = "; ")
      ), call. = FALSE)
    }
  } else {
    for (i in names(list_train_vali_Data)) {
      common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
    }
    dropped <- setdiff(candidate_genes, common_feature[-c(1:2)])
    if (length(dropped) > 0L) {
      warning(paste0(
        "feature_alignment='intersection' dropped candidate features absent from at least one cohort: ",
        paste(dropped, collapse = ", ")
      ), call. = FALSE)
    }
  }

  ##### parameters check #####


  validate_category_ml_inputs(
    train_data = train_data,
    list_train_vali_Data = list_train_vali_Data,
    candidate_genes = candidate_genes,
    common_feature = common_feature
  )
  methods <- resolve_category_methods(methods)
  if (length(cores_for_parallel) != 1L || is.na(cores_for_parallel) ||
      cores_for_parallel < 1L) {
    stop("cores_for_parallel must be a positive scalar integer", call. = FALSE)
  }
  cores_for_parallel <- as.integer(cores_for_parallel)

  if (
    #identical(colnames(train_data)[1:2], c("ID", "Var")) &
    #identical(common_feature[1:2], c("ID", "Var")) &
    #unique(train_data$Var) %in% c("Y", "N") &
    #length(candidate_genes) > 1 &
	    all(is.element(methods, c("nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass")))


  ) {
    ####### data preparation ######





	    category_prepped <- fit_category_preprocess_recipe(
	      train_data,
	      common_feature = common_feature,
	      label = "Training data",
	      positive_class = positive_class
	    )
	    train_data <- category_prepped$data
	    preprocess_recipe <- category_prepped$recipe
	    list_train_vali_Data <- lapply(names(list_train_vali_Data), function(nm) {
	      apply_category_preprocess_recipe(
	        list_train_vali_Data[[nm]],
	        recipe = preprocess_recipe,
	        label = paste0("Dataset '", nm, "'"),
	        common_feature = common_feature
	      )
	    })
	    names(list_train_vali_Data) <- data_names

	    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
	    pre_var <- preprocess_recipe$feature_names



    print(paste0("There existing ", length(candidate_genes), " genes in candidate genes"))


    print("Intersetion of the candidate genes and the colnames of the provided data")

    print(paste0("There existing ", length(pre_var), " genes in candidate genes, colnames of training data, colnames of validation data"))




    # parallel processing


	  cl <- parallel::makePSOCKcluster(cores_for_parallel)
	    on.exit({
	      if (!is.null(cl)) {
	        parallel::stopCluster(cl)
	      }
	      foreach::registerDoSEQ()
	    }, add = TRUE)
	    doParallel::registerDoParallel(cl)

    set.seed(seed)
    res.model <- model.Dev(
      training = train_data,
      method = methods,
      sig = pre_var
    )
	    parallel::stopCluster(cl)
	    cl <- NULL
	    foreach::registerDoSEQ()



    ml.auc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.auc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.auc) <- names(list_train_vali_Data)

    ml.roc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.roc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.roc) <- names(list_train_vali_Data)

    res <- list()
	    res[["model"]] <- res.model
	    res[["auc"]] <- ml.auc
	    res[["roc"]] <- ml.roc
	    res[["sig.gene"]] <- pre_var
	    res[["Preprocess.recipe"]] <- preprocess_recipe
	    res[["positive_class"]] <- positive_class

	    return(res)
  } else {
    print("Please provide the correct parameters")
  }
}
