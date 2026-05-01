#' Screening out the core variables for the prognosis with the machine learning algorithms
#' 
#' A function can be used for screening out the core features from the given candidate genes with seven machine learning algorithms.
#' 
#' @param InputMatrix A gene expression dataframe after log2(x+1) scaled. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead. 
#' @param candidate_genes The input genes, that you want to screen out from, for identifying the core features.
#' @param mode  We provide three modes including 'all', 'single', and 'all_without_SVM'. The 'all' mode means using all seven methods for selecting. The 'single' mode means using only one method for running. Since SVM takes so much time, we're singling him out. The 'all_without_SVM' mode means the other six methods used for selecting.
#' @param seed  The seed. You can set it as any number. For example, 5201314.
#' @param single_ml The one method from the seven methods including "RSF", "Enet", "Xgboost", "SVM-REF", "Lasso", "CoxBoost", "StepCox".
#' @param nodesize The node size parameter for the RSF method. The default is 5. You can try another positive integer. For example, 10,15,20, etc. 
#'
#' @return A data frame including the methods and the core genes screened by the corresponding algorithm.
#' @export
#'
#' @examples
ML.Corefeature.Prog.Screen <- function(InputMatrix, ### ???ID,???OS.time, (day), ??? OS, (0/1), ?????????????????log2(x+1)
                                       candidate_genes,
                                       mode = NULL, # all, single,all_without_SVM
                                       seed = NULL,
                                       single_ml = NULL, # c("RSF", "Enet", "Xgboost","SVM-REF","Lasso","CoxBoost","StepCox")
                                       nodesize = 5) {
  ### laoding the function ####

  if (T) {
    svmRFE.wrap <- function(test.fold, X, ...) {
      # Wrapper to run svmRFE function while omitting a given test fold
      train.data <- X[-test.fold, ]
      test.data <- X[test.fold, ]

      # Rank the features
      features.ranked <- svmRFE(train.data, ...)

      return(list(feature.ids = features.ranked, train.data.ids = row.names(train.data), test.data.ids = row.names(test.data)))
    }

    svmRFE <- function(X, k = 1, halve.above = 5000) {
      # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
      n <- ncol(X) - 1

      # Scale data up front so it doesn't have to be redone each pass
      cat("Scaling data...")
      X[, -1] <- scale(X[, -1])
      cat("Done!\n")
      flush.console()

      pb <- txtProgressBar(1, n, 1, style = 3)

      i.surviving <- 1:n
      i.ranked <- n
      ranked.list <- vector(length = n)

      # Recurse through all the features
      while (length(i.surviving) > 0) {
        if (k > 1) {
          # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)
          folds <- rep(1:k, len = nrow(X))[sample(nrow(X))]
          folds <- lapply(1:k, function(x) which(folds == x))

          # Obtain weights for each training set
          w <- lapply(folds, getWeights, X[, c(1, 1 + i.surviving)])
          w <- do.call(rbind, w)

          # Normalize each weights vector
          w <- t(apply(w, 1, function(x) x / sqrt(sum(x^2))))

          # Compute ranking criteria
          v <- w * w
          vbar <- apply(v, 2, mean)
          vsd <- apply(v, 2, sd)
          c <- vbar / vsd
        } else {
          # Only do 1 pass (i.e. regular SVM-RFE)
          w <- getWeights(NULL, X[, c(1, 1 + i.surviving)])
          c <- w * w
        }

        # Rank the features
        ranking <- sort(c, index.return = T)$ix
        if (length(i.surviving) == 1) {
          ranking <- 1
        }

        if (length(i.surviving) > halve.above) {
          # Cut features in half until less than halve.above
          nfeat <- length(i.surviving)
          ncut <- round(nfeat / 2)
          n <- nfeat - ncut

          cat("Features halved from", nfeat, "to", n, "\n")
          flush.console()

          pb <- txtProgressBar(1, n, 1, style = 3)
        } else {
          ncut <- 1
        }

        # Update feature list
        ranked.list[i.ranked:(i.ranked - ncut + 1)] <- i.surviving[ranking[1:ncut]]
        i.ranked <- i.ranked - ncut
        i.surviving <- i.surviving[-ranking[1:ncut]]

        setTxtProgressBar(pb, n - length(i.surviving))
        flush.console()
      }

      close(pb)

      return(ranked.list)
    }

    getWeights <- function(test.fold, X) {
      # Fit a linear SVM model and obtain feature weights
      train.data <- X
      if (!is.null(test.fold)) train.data <- X[-test.fold, ]

      svmModel <- e1071::svm(train.data[, -1], train.data[, 1],
        cost = 10, cachesize = 500,
        scale = F, type = "C-classification", kernel = "linear"
      )

      t(svmModel$coefs) %*% svmModel$SV
    }

    survival_svm_screen <- function(est_dd, seed, max_features = 300) {
      features <- colnames(est_dd)[-c(1, 2)]
      if (length(features) <= 1L) {
        return(features)
      }
      fit <- train_survivalsvm(est_dd, seed)
      base_rs <- predict_survivalsvm(fit, est_dd)
      base_cindex <- calculate_cindex(base_rs, est_dd)
      if (is.na(base_cindex)) {
        base_cindex <- 0.5
      }

      importance <- vapply(seq_along(features), function(i) {
        feature <- features[[i]]
        permuted <- est_dd
        set.seed(seed + i)
        permuted[[feature]] <- sample(permuted[[feature]])
        rs <- tryCatch(predict_survivalsvm(fit, permuted), error = function(e) rep(NA_real_, nrow(permuted)))
        cindex <- calculate_cindex(rs, permuted)
        if (is.na(cindex)) {
          return(NA_real_)
        }
        base_cindex - cindex
      }, numeric(1))

      names(importance) <- features
      ordered <- names(sort(importance, decreasing = TRUE, na.last = NA))
      selected <- ordered[importance[ordered] > 0]
      if (length(selected) < 2L) {
        selected <- ordered
      }
      head(selected, min(max_features, length(selected)))
    }

	    survival_xgboost_screen <- function(est_dd, seed) {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop(
          "ML.Corefeature.Prog.Screen single_ml='Xgboost' or Xgboost-containing modes require the optional xgboost package; install xgboost or choose a non-Xgboost mode.",
          call. = FALSE
        )
      }
      features <- colnames(est_dd)[-c(1, 2)]
      train <- as.data.frame(est_dd[, features, drop = FALSE])
      train_matrix <- Matrix::sparse.model.matrix(~ . - 1, data = train)
      # xgboost survival:cox treats negative labels as right-censored times.
      train_label <- as.numeric(est_dd$OS.time)
      train_label[as.numeric(est_dd$OS) == 0] <- -train_label[as.numeric(est_dd$OS) == 0]
      dtrain <- xgboost::xgb.DMatrix(data = train_matrix, label = train_label)
      xgb <- xgboost::xgboost(
        data = dtrain,
        max_depth = 3,
        eta = 0.05,
        objective = "survival:cox",
        nrounds = 100,
        verbose = 0
      )
      importance <- xgboost::xgb.importance(colnames(train_matrix), model = xgb)
      if (nrow(importance) == 0L) {
        return(character())
      }
	      importance$rel.imp <- importance$Gain / max(importance$Gain, na.rm = TRUE)
	      importance[importance$rel.imp >= 0.05, "Feature"]
	    }

	    append_screen_result <- function(selected.feature, method, rid) {
	      result <- data.frame(
	        method = c(rep(method, length(rid))),
	        selected.fea = rid
	      )
	      rbind(selected.feature, result)
	    }

	    screen_rsf_vars <- function(est_dd, rf_nodesize, seed) {
	      fit <- train_rsf(est_dd, rf_nodesize = rf_nodesize, seed = seed)
	      get_rsf_selected_vars(fit, seed = seed)
	    }

	    screen_coxboost_vars <- function(est_dd, seed) {
	      fit <- train_coxboost(est_dd, seed = seed)
	      get_coxboost_selected_vars(fit)
	    }

	    screen_stepcox_vars <- function(est_dd, direction) {
	      fit <- train_stepcox(est_dd, direction = direction)
	      get_stepcox_selected_vars(fit)
	    }

	    WriteFeatures <- function(results, input, save = T, file = "features_ranked.txt") {
      # Compile feature rankings across multiple folds
      featureID <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$ix
      avg.rank <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$x
      feature.name <- colnames(input[, -1])[featureID]
      features.ranked <- data.frame(FeatureName = feature.name, FeatureID = featureID, AvgRank = avg.rank)
      if (save == T) {
        write.table(features.ranked, file = file, quote = F, row.names = F)
      } else {
        features.ranked
      }
    }

    FeatSweep.wrap <- function(i, results, input) {
      # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
      svm.list <- lapply(results, function(x) {
        e1071::tune(e1071::svm,
          train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
          train.y = input[x$train.data.ids, 1],
          validation.x = input[x$test.data.ids, 1 + x$feature.ids[1:i]],
          validation.y = input[x$test.data.ids, 1],
          # Optimize SVM hyperparamters
          ranges = e1071::tune(e1071::svm,
            train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
            train.y = input[x$train.data.ids, 1],
            ranges  = list(gamma = 2^(-12:0), cost = 2^(-6:6))
          )$best.par,
          tunecontrol = e1071::tune.control(sampling = "fix")
        )$perf
      })

      error <- mean(sapply(svm.list, function(x) x$error))
      return(list(svm.list = svm.list, error = error))
    }


    PlotErrors <- function(errors, errors2 = NULL, no.info = 0.5,
                           ylim = range(c(errors, errors2), na.rm = T),
                           xlab = "Number of Features", ylab = "5 x CV Error") {
      # Makes a plot of average generalization error vs. number of top features
      AddLine <- function(x, col = "dodgerblue") {
        lines(which(!is.na(errors)), na.omit(x), col = col, lwd = 3)
        points(which.min(x), min(x, na.rm = T), col = "firebrick3")
        text(which.min(x), min(x, na.rm = T), paste(
          which.min(x), "-",
          format(min(x, na.rm = T), dig = 3)
        ), pos = 2, col = "red", cex = 1.15)
      }

      plot(errors, type = "n", ylim = ylim, xlab = xlab, ylab = ylab)
      AddLine(errors)
      if (!is.null(errors2)) AddLine(errors2, "gray30")
      abline(h = no.info, lty = 2)
    }


    Plotaccuracy <- function(errors, errors2 = NULL, no.info = 0.5,
                             ylim = range(c(errors, errors2), na.rm = T),
                             xlab = "Number of Features", ylab = "5 x CV Accuracy") {
      # Makes a plot of average generalization error vs. number of top features
      AddLine <- function(x, col = "dodgerblue") {
        lines(which(!is.na(errors)), na.omit(x), col = col, lwd = 3)
        points(which.max(x), max(x, na.rm = T), col = "firebrick3")
        text(which.max(x), max(x, na.rm = T), paste(
          which.max(x), "-",
          format(max(x, na.rm = T), dig = 3)
        ), pos = 2, col = "red", cex = 1.15)
      }

      plot(errors, type = "n", ylim = ylim, xlab = xlab, ylab = ylab)
      AddLine(errors)
      if (!is.null(errors2)) AddLine(errors2, "gray30")
      abline(h = no.info, lty = 2)
    }


    prepare_core_screening_input <- function(gene_list, inputSet, label) {
      print("Starting the data preprocess")
      print("Rejecting a null value")

      gene_list <- normalize_ml_feature_names(gene_list)
      inputSet <- normalize_ml_data_columns(inputSet, label)

      print("Gets the intersection of genelist and expression profile")
      common_feature <- intersect(c("ID", "OS.time", "OS", gene_list), colnames(inputSet))
      missing_required <- setdiff(c("ID", "OS.time", "OS"), common_feature)
      if (length(missing_required) > 0L) {
        stop(paste0(label, " is missing required columns: ",
                    paste(missing_required, collapse = ", ")), call. = FALSE)
      }
      if (length(common_feature) <= 3L) {
        stop(paste0(label, " has no candidate genes present after name normalization"),
             call. = FALSE)
      }

      dropped <- setdiff(gene_list, common_feature[-c(1:3)])
      if (length(dropped) > 0L) {
        warning(paste0(label, " dropped candidate genes absent from the expression matrix: ",
                       paste(dropped, collapse = ", ")), call. = FALSE)
      }

      print("Processing the input representation matrix")
      prepped <- fit_survival_preprocess_recipe(inputSet, common_feature, label)$data
      print("Data preprocessing completed")
      prepped
    }

    screening_feature_names <- function(inputSet) {
      setdiff(colnames(inputSet), c("ID", "OS.time", "OS"))
    }

    SigKMcox <- function(gene_list,
                         inputSet,
                         KM_pcutoff # KM?????
    ) {
      inputSet <- prepare_core_screening_input(gene_list, inputSet, "ML.Corefeature.Prog.Screen KM prefilter")
      features <- screening_feature_names(inputSet)

      display.progress <- function(index, totalN, breakN = 20) {
        if (totalN > 0L && index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }

      print("Stating the KM selection")
      kmoutput <- data.frame()

      for (i in seq_along(features)) {
        display.progress(index = i, totalN = length(features), breakN = 20)
        g <- features[[i]]
        tmp <- inputSet[, c("OS.time", "OS", g), drop = FALSE]
        if (length(unique(tmp[[g]])) < 2L) {
          next
        }
        group <- ifelse(tmp[[g]] > stats::median(tmp[[g]]), "High", "Low")
        if (length(unique(group)) < 2L) {
          next
        }
        tmp$group <- factor(group, levels = c("Low", "High"))
        p.val <- tryCatch({
          fitd <- survival::survdiff(
            survival::Surv(OS.time, OS) ~ group,
            data = tmp,
            na.action = stats::na.exclude
          )
          stats::pchisq(fitd$chisq, length(fitd$n) - 1, lower.tail = FALSE)
        }, error = function(e) NA_real_)
        kmoutput <- rbind(kmoutput, data.frame(
          gene = g,
          pvalue = p.val,
          stringsAsFactors = FALSE
        ))
      }

      print("Finished the KM selection")
      if (nrow(kmoutput) == 0L) {
        return(character())
      }
      kmoutput[which(!is.na(kmoutput$pvalue) & kmoutput$pvalue < KM_pcutoff), "gene"]
    }


    SigUnicox <- function(gene_list,
                          inputSet,
                          unicox_pcutoff # ??????????
    ) {
      inputSet <- prepare_core_screening_input(gene_list, inputSet, "ML.Corefeature.Prog.Screen unicox prefilter")
      features <- screening_feature_names(inputSet)

      display.progress <- function(index, totalN, breakN = 20) {
        if (totalN > 0L && index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }

      print("Stating the univariable cox regression")
      unicox <- data.frame()
      for (i in seq_along(features)) {
        display.progress(index = i, totalN = length(features))
        gene <- features[[i]]
        tmp <- data.frame(
          expr = as.numeric(inputSet[[gene]]),
          futime = inputSet$OS.time,
          fustat = inputSet$OS,
          stringsAsFactors = FALSE
        )
        if (length(unique(tmp$expr)) < 2L) {
          next
        }
        coxSummary <- tryCatch({
          cox <- survival::coxph(survival::Surv(futime, fustat) ~ expr, data = tmp)
          summary(cox)
        }, error = function(e) NULL)
        if (is.null(coxSummary) || nrow(coxSummary$coefficients) == 0L) {
          next
        }
        unicox <- rbind.data.frame(unicox,
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

      print("Finished the univariable cox regression")
      if (nrow(unicox) == 0L) {
        return(character())
      }
      unicox[which(!is.na(unicox$pvalue) & unicox$pvalue < unicox_pcutoff), "gene"]
    }
  }

  valid_core_screen_methods <- c("RSF", "Enet", "Xgboost", "SVM-REF",
                                 "Lasso", "CoxBoost", "StepCox")
  if (identical(single_ml, "Boruta")) {
    stop(
      "Boruta core-feature screening has been removed because it duplicated RSF ",
      "after the survival-aware fix. Use single_ml = 'RSF' for survival-aware ",
      "random-forest feature screening.",
      call. = FALSE
    )
  }
  if (identical(mode, "single") &&
      (is.null(single_ml) || length(single_ml) != 1L ||
       is.na(single_ml) || !single_ml %in% valid_core_screen_methods)) {
    stop(paste0(
      "single_ml must be one of: ",
      paste(valid_core_screen_methods, collapse = ", ")
    ), call. = FALSE)
  }


  needs_xgboost <- identical(mode, "all") ||
    identical(mode, "all_without_SVM") ||
    identical(single_ml, "Xgboost")
  if (needs_xgboost) {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop(
        "ML.Corefeature.Prog.Screen single_ml='Xgboost' or Xgboost-containing modes require the optional xgboost package; install xgboost or choose a non-Xgboost mode.",
        call. = FALSE
      )
    }
  }

  InputMatrix_pre <- InputMatrix
  colnames(InputMatrix_pre) <- gsub("-", ".", colnames(InputMatrix_pre))
  genelist.1 <- SigUnicox(gene_list = candidate_genes, inputSet = InputMatrix_pre, unicox_pcutoff = 0.05)
  if (length(genelist.1) < 1L) {
    stop(
      "No genes passed the univariate Cox prefilter in ML.Corefeature.Prog.Screen; ",
      "relax the cutoff or provide candidate genes with stronger survival signal.",
      call. = FALSE
    )
  }

  genelist.2 <- SigKMcox(gene_list = genelist.1, inputSet = InputMatrix_pre, KM_pcutoff = 0.05)
  if (length(genelist.2) < 1L) {
    stop(
      "No genes passed the KM prefilter in ML.Corefeature.Prog.Screen; ",
      "relax the cutoff or provide candidate genes with stronger survival signal.",
      call. = FALSE
    )
  }

  candidate_genes <- genelist.2

  print("----- finish the preprocess of the unicox and km analysis-----")

  ##### setting the pamameters ######

  rf_nodesize <- nodesize
  seed <- seed
  iter.times <- 1000
  # Checking data feasibility
  message("--- check data feasibility ---")

  # Replace '-' in column names with '.'

  candidate_genes <- normalize_ml_feature_names(candidate_genes)
  colnames(InputMatrix) <- normalize_ml_feature_names(colnames(InputMatrix))


  # Matching candidate genes to genes in each cohort
  common_feature <- c("ID", "OS.time", "OS", candidate_genes)
  common_feature <- intersect(common_feature, colnames(InputMatrix))

  message(paste0("---the number of the raw candidate genes is ", length(candidate_genes), " ---"))
  message(paste0("---the number of the common feature is ", length(common_feature) - 3, " ---"))

  ######### the main of the function ##########

  if (!is.na(rf_nodesize) &
    !is.na(seed) &
    mode %in% c("all", "single", "all_without_SVM") &
    identical(c("ID", "OS.time", "OS"), colnames(InputMatrix)[1:3]) &
    length(candidate_genes) > 0 &
    identical(c("ID", "OS.time", "OS"), common_feature[1:3]) &
    length(common_feature) > 3) {
    message("--- Data preprocessing ---")
    # Data preprocessing

    # Reuse the survival ML preprocessing contract so screening does not diverge
    # from model development on NA handling, event coding, or ID validation.
    InputMatrix <- fit_survival_preprocess_recipe(
      InputMatrix,
      common_feature,
      label = "ML.Corefeature.Prog.Screen input"
    )$data

    est_dd <- as.data.frame(InputMatrix)[, common_feature[-1]]
    pre_var <- common_feature[-c(1:3)]
    selected.feature <- data.frame()


    if (mode == "all") {
      ### 1. Repeated Lasso  #############
      message("--- 1.Repeated lasso ---")
      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("1000 time lasso penalty")
      # 1000 time lasso penalty
      lasso_fea_list <- list()
      list.of.seed <- 1:iter.times

      print("This step will probably take several hours")

      lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
        set.seed(list.of.seed[x])
        cvfit <- glmnet::cv.glmnet(
          x = x1,
          y = x2,
          nfolds = 10, # 10-fold????????lambda
          alpha = 1, # alpha = 1 ??? lasso
          family = "cox", # ??cox??
          maxit = 1000
        )

        # optimal lambda
        fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
        if (is.element("(Intercept)", fea)) {
          lasso_fea <- sort(fea[-1]) # ????????
        } else {
          lasso_fea <- sort(fea)
        }
        return(lasso_fea)
      })

      # ???????????
      lasso_res <- NULL
      for (i in 1:iter.times) {
        lasso_res <- rbind.data.frame(lasso_res,
          data.frame(
            iteration = i,
            n.gene = length(lasso_fea_list[[i]]),
            genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
            stringsAsFactors = F
          ),
          stringsAsFactors = F
        )
      }




      genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # ???????????
      freq.cutoff <- iter.times * 0.05
      genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????. 95%


      result <- data.frame(
        method = c(rep("Lasso", length(genes))),
        selected.fea = genes
      )

      selected.feature <- rbind(selected.feature, result)



      ##### 2.Enet ###########
      message("--- 2.Enet  ---")


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("This step will probably take several hours")

      for (alpha in seq(0.1, 0.9, 0.1)) {
        message(paste0("--- 2.Enet ", alpha, "---"))

        set.seed(seed)
        # 1000 time
        fea_list <- list()
        list.of.seed <- 1:iter.times
        fea_list <- pblapply(list.of.seed, function(x) { # ????2?
          set.seed(list.of.seed[x])
          cvfit <- glmnet::cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold????????lambda
            alpha = alpha,
            family = "cox", # ??cox??
            maxit = 1000
          )
          # ????lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # ????????
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })


        genes <- sort(table(unlist(fea_list)), decreasing = T) # ???????????
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????

        result <- data.frame(
          method = c(rep(paste0("Enet", "[?=", alpha, "]"), length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)
      }

      ##### 3.SVM-REF ###########
      message("--- 3.SVM-REF  ---")
      message("Using survival-SVM permutation importance with OS.time/OS")
      fea <- survival_svm_screen(est_dd, seed)

      result <- data.frame(
        method = c(rep("SVM-REF", length(fea))),
        selected.fea = fea
      )

      selected.feature <- rbind(selected.feature, result)

      ##### 4.Xgboost ###########

      message("--- 4.Xgboost  ---")
      xgboost.finalVars <- survival_xgboost_screen(est_dd, seed)

      result <- data.frame(
        method = c(rep("Xgboost", length(xgboost.finalVars))),
        selected.fea = xgboost.finalVars
      )

      selected.feature <- rbind(selected.feature, result)


	      ##### 5.RSF ###########
	      message("--- 5.RSF  ---")

	      rid <- screen_rsf_vars(est_dd, rf_nodesize, seed)
	      selected.feature <- append_screen_result(selected.feature, "RSF", rid)


	      ##### 6.CoxBoost ###########
	      message("--- 6.CoxBoost  ---")

	      rid <- screen_coxboost_vars(est_dd, seed)
	      selected.feature <- append_screen_result(selected.feature, "CoxBoost", rid)

      ##### 7.StepCox ###########
	      message("--- 7.StepCox ---")

	      for (direction in c("both", "backward", "forward")) {
	        rid <- screen_stepcox_vars(est_dd, direction)
	        selected.feature <- append_screen_result(
	          selected.feature,
	          paste0("StepCox", "+", direction),
	          rid
	        )
	      }

      return(selected.feature)
    } else if (mode == "single") {
      if (single_ml == "Lasso") {
        ### 1. Repeated Lasso  #############
        message("--- 1.Repeated lasso ---")
        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        print("1000 time lasso penalty")
        # 1000 time lasso penalty
        lasso_fea_list <- list()
        list.of.seed <- 1:iter.times

        print("This step will probably take several hours")

        lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
          set.seed(list.of.seed[x])
          cvfit <- glmnet::cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold????????lambda
            alpha = 1, # alpha = 1 ??? lasso
            family = "cox", # ??cox??
            maxit = 1000
          )

          # optimal lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # ????????
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })

        # ???????????
        lasso_res <- NULL
        for (i in 1:iter.times) {
          lasso_res <- rbind.data.frame(lasso_res,
            data.frame(
              iteration = i,
              n.gene = length(lasso_fea_list[[i]]),
              genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
              stringsAsFactors = F
            ),
            stringsAsFactors = F
          )
        }




        genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # ???????????
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????. 95%


        result <- data.frame(
          method = c(rep("Lasso", length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "Enet") {
        ##### 2.Enet ###########
        message("--- 2.Enet  ---")


        x1 <- as.matrix(est_dd[, pre_var])
        x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
        print("This step will probably take several hours")

        for (alpha in seq(0.1, 0.9, 0.1)) {
          message(paste0("--- 2.Enet ", alpha, "---"))

          set.seed(seed)
          # 1000 time
          fea_list <- list()
          list.of.seed <- 1:iter.times
          fea_list <- pblapply(list.of.seed, function(x) { # ????2?
            set.seed(list.of.seed[x])
            cvfit <- glmnet::cv.glmnet(
              x = x1,
              y = x2,
              nfolds = 10, # 10-fold????????lambda
              alpha = alpha,
              family = "cox", # ??cox??
              maxit = 1000
            )
            # ????lambda
            fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
            if (is.element("(Intercept)", fea)) {
              lasso_fea <- sort(fea[-1]) # ????????
            } else {
              lasso_fea <- sort(fea)
            }
            return(lasso_fea)
          })


          genes <- sort(table(unlist(fea_list)), decreasing = T) # ???????????
          freq.cutoff <- iter.times * 0.05
          genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????

          result <- data.frame(
            method = c(rep(paste0("Enet", "[?=", alpha, "]"), length(genes))),
            selected.fea = genes
          )

          selected.feature <- rbind(selected.feature, result)
        }

        return(selected.feature)
      } else if (single_ml == "SVM-REF") {
        ##### 3.SVM-REF ###########
        message("--- 3.SVM-REF  ---")
        message("Using survival-SVM permutation importance with OS.time/OS")
        fea <- survival_svm_screen(est_dd, seed)

        result <- data.frame(
          method = c(rep("SVM-REF", length(fea))),
          selected.fea = fea
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
      } else if (single_ml == "Xgboost") {
        ##### 4.Xgboost ###########

        message("--- 4.Xgboost  ---")
        xgboost.finalVars <- survival_xgboost_screen(est_dd, seed)

        result <- data.frame(
          method = c(rep("Xgboost", length(xgboost.finalVars))),
          selected.fea = xgboost.finalVars
        )

        selected.feature <- rbind(selected.feature, result)

        return(selected.feature)
	      } else if (single_ml == "RSF") {
	        ##### 5.RSF ###########
	        message("--- 5.RSF  ---")

	        rid <- screen_rsf_vars(est_dd, rf_nodesize, seed)
	        selected.feature <- append_screen_result(selected.feature, "RSF", rid)

	        return(selected.feature)
	      } else if (single_ml == "CoxBoost") {
	        ##### 6.CoxBoost ###########
	        message("--- 6.CoxBoost  ---")

	        rid <- screen_coxboost_vars(est_dd, seed)
	        selected.feature <- append_screen_result(selected.feature, "CoxBoost", rid)
	        return(selected.feature)
	      } else if (single_ml == "StepCox") {
	        ##### 7.StepCox ###########
	        message("--- 7.StepCox ---")

	        for (direction in c("both", "backward", "forward")) {
	          rid <- screen_stepcox_vars(est_dd, direction)
	          selected.feature <- append_screen_result(
	            selected.feature,
	            paste0("StepCox", "+", direction),
	            rid
	          )
	        }

        return(selected.feature)
      } else {
        warning("The parameter of the single ML is out of the bound")
      }
    } else if (mode == "all_without_SVM") {
      ### 1. Repeated Lasso  #############
      message("--- 1.Repeated lasso ---")
      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("1000 time lasso penalty")
      # 1000 time lasso penalty
      lasso_fea_list <- list()
      list.of.seed <- 1:iter.times

      print("This step will probably take several hours")

      lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
        set.seed(list.of.seed[x])
        cvfit <- glmnet::cv.glmnet(
          x = x1,
          y = x2,
          nfolds = 10, # 10-fold????????lambda
          alpha = 1, # alpha = 1 ??? lasso
          family = "cox", # ??cox??
          maxit = 1000
        )

        # optimal lambda
        fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
        if (is.element("(Intercept)", fea)) {
          lasso_fea <- sort(fea[-1]) # ????????
        } else {
          lasso_fea <- sort(fea)
        }
        return(lasso_fea)
      })

      # ???????????
      lasso_res <- NULL
      for (i in 1:iter.times) {
        lasso_res <- rbind.data.frame(lasso_res,
          data.frame(
            iteration = i,
            n.gene = length(lasso_fea_list[[i]]),
            genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
            stringsAsFactors = F
          ),
          stringsAsFactors = F
        )
      }




      genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # ???????????
      freq.cutoff <- iter.times * 0.05
      genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????. 95%


      result <- data.frame(
        method = c(rep("Lasso", length(genes))),
        selected.fea = genes
      )

      selected.feature <- rbind(selected.feature, result)



      ##### 2.Enet ###########
      message("--- 2.Enet  ---")


      x1 <- as.matrix(est_dd[, pre_var])
      x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
      print("This step will probably take several hours")

      for (alpha in seq(0.1, 0.9, 0.1)) {
        message(paste0("--- 2.Enet ", alpha, "---"))

        set.seed(seed)
        # 1000 time
        fea_list <- list()
        list.of.seed <- 1:iter.times
        fea_list <- pblapply(list.of.seed, function(x) { # ????2?
          set.seed(list.of.seed[x])
          cvfit <- glmnet::cv.glmnet(
            x = x1,
            y = x2,
            nfolds = 10, # 10-fold????????lambda
            alpha = alpha,
            family = "cox", # ??cox??
            maxit = 1000
          )
          # ????lambda
          fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
          if (is.element("(Intercept)", fea)) {
            lasso_fea <- sort(fea[-1]) # ????????
          } else {
            lasso_fea <- sort(fea)
          }
          return(lasso_fea)
        })


        genes <- sort(table(unlist(fea_list)), decreasing = T) # ???????????
        freq.cutoff <- iter.times * 0.05
        genes <- names(genes[genes > freq.cutoff]) # ??????????50?????????lasso?????

        result <- data.frame(
          method = c(rep(paste0("Enet", "[?=", alpha, "]"), length(genes))),
          selected.fea = genes
        )

        selected.feature <- rbind(selected.feature, result)
      }

      ##### 3.SVM-REF ###########
      # message('--- 3.SVM-REF  ---')
      # print('This step will probably take several hours')
      #
      #
      # input <- est_dd[,-1]
      #
      #
      #
      # # 10CV (k-fold crossValidation?
      # svmRFE(input, k = 10, halve.above = 100) #??????????
      # nfold = 10
      # nrows = nrow(input)
      # folds = rep(1:nfold, len=nrows)[sample(nrows)]
      # folds = lapply(1:nfold, function(x) which(folds == x))
      # results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) #????
      # top.features = WriteFeatures(results, input, save=F) #??????
      # n.features = nrow(top.features)
      # if(n.features > 300){
      #   n.svm = 300
      # }else {
      #   n.svm =n.features
      # }
      #
      # featsweep = base::lapply(1:n.svm, FeatSweep.wrap, results, input)
      #
      # no.info = min(prop.table(table(input[,1])))
      # errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
      # fea <- top.features[1:which.min(errors), "FeatureName"]
      #
      # result <-  data.frame(method = c(rep("SVM-REF", length(fea))),
      #                       selected.fea= fea)
      #
      # selected.feature <- rbind(selected.feature,result)

      ##### 3.Xgboost ###########

      message("--- 3.Xgboost  ---")
      xgboost.finalVars <- survival_xgboost_screen(est_dd, seed)

      result <- data.frame(
        method = c(rep("Xgboost", length(xgboost.finalVars))),
        selected.fea = xgboost.finalVars
      )

      selected.feature <- rbind(selected.feature, result)


	      ##### 4.RSF ###########
	      message("--- 4.RSF  ---")

	      rid <- screen_rsf_vars(est_dd, rf_nodesize, seed)
	      selected.feature <- append_screen_result(selected.feature, "RSF", rid)


	      ##### 5.CoxBoost ###########
	      message("--- 5.CoxBoost  ---")

	      rid <- screen_coxboost_vars(est_dd, seed)
	      selected.feature <- append_screen_result(selected.feature, "CoxBoost", rid)

      ##### 6.StepCox ###########
	      message("--- 6.StepCox ---")

	      for (direction in c("both", "backward", "forward")) {
	        rid <- screen_stepcox_vars(est_dd, direction)
	        selected.feature <- append_screen_result(
	          selected.feature,
	          paste0("StepCox", "+", direction),
	          rid
	        )
	      }

      return(selected.feature)
    }
  } else {
    print("Please set the full parameters")
  }
}
