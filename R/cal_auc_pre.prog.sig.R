#' Calculate AUC scores of the previous signarues in all data
#'
#' @param use_your_own_collected_sig Whether to use your own collected signatures. T or F.
#' @param collected_sig_table If use_your_own_collected_sig set as T, you should provide a data frame containing the information of the signatures. The column names of the data frame are "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol". For example, 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'. The 'model' consists of the first name of the first author and the PMID. The 'PMID' is from the paper. The 'Cancer' uses abbreviations in the format of the TCGA. 'Author' is the first name of the first author of the paper. 'Coef' is the coefficient of the variable in the signature. 'symbol' is the variable in the signature. If use_your_own_collected_sig is F, you don't need to provide this data.
#' @param type.sig If the use_your_own_collected_sig is F, here we collected some signatures of the LGG, GBM, and Glioma. You can choose c('Glioma', 'LGG', 'GBM'), c('Glioma'),c('GBM'), c('Glioma', 'LGG'), and some other combination for the signatures you want.
#' @param list_input_data A list of the cohorts. Column names are 'ID', 'OS.time', 'OS', and the other variables. OS.time (Day). OS(1/0).
#' @param AUC_time Time horizon in years, for example 1, 2, or 3. Use a horizon shorter than the minimum maximum survival time across cohorts.
#' @param auc_cal_method 'KM', 'NNE'. The default is 'KM'.
#'
#' @return A list of the AUC results of each previous signature in each cohort you provide.
#' @export
#'
cal_auc_pre.prog.sig <- function(use_your_own_collected_sig, # 是否使用您自己收集的signature， T or F
                                 collected_sig_table, # 列名分别为
                                 # "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol"
                                 # 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'
                                 type.sig, ### prognostic signature 的类型，c('Glioma','LGG','GBM')， c('Glioma'),c('Glioma','LGG')
                                 list_input_data, # list of the cohorts(ID,OS.time, OS····)经过了log2（x+1）转化
                                 AUC_time = NULL, ### 时间 年份
                                 auc_cal_method = "KM" # KM, NNE 默认为KM
) {
  if (use_your_own_collected_sig) {
    sig.input <- collected_sig_table
  } else {
    pre_prog_sig_data <- get("pre.prog.sig", envir = asNamespace("iklSurvML"))

    if (all(type.sig %in% names(pre_prog_sig_data))) {
      if (length(type.sig) == 1) {
        sig.input <- pre_prog_sig_data[[type.sig[1]]]
      } else {
        sig.input <- pre_prog_sig_data[[type.sig[1]]]

        for (i in 2:length(type.sig)) {
          sig.input <- rbind(sig.input, pre_prog_sig_data[[type.sig[i]]])
        }
      }
    } else {
      stop("please provide correct type.sig")
    }
  }





  sig.input$Coef <- as.numeric(sig.input$Coef)
  sig.input$symbol <- gsub("-", ".", sig.input$symbol)



  # Replace '-' in column names with '.'
  list_input_data <- lapply(list_input_data, function(x) {
    colnames(x) <- gsub("-", ".", colnames(x))
    return(x)
  })

  common_feature <- c("ID", "OS.time", "OS", unique(sig.input$symbol))

  # for (i in names(list_input_data)) {
  #   common_feature = intersect(common_feature, colnames(list_input_data[[i]]))
  # }

  returnIDtoRS <- function(rs.table.list, rawtableID) {
    for (i in names(rs.table.list)) {
      rs.table.list[[i]]$ID <- rawtableID[[i]]$ID
      rs.table.list[[i]] <- rs.table.list[[i]] %>% dplyr::select("ID", everything())
    }

    return(rs.table.list)
  }


  returnRStoROC <- function(rs.table.list, AUC_time) {
    lapply(rs.table.list, function(x) {
      if (length(unique(x$RS)) <= 1L) {
        return(data.frame(
          TP = rep(0, nrow(x)),
          FP = rep(1, nrow(x)),
          AUC = rep(0.5, nrow(x)),
          HR = rep(1, nrow(x)),
          marker_direction = "higher_is_worse"
        ))
      }
      calculate_survival_roc_from_risk(
        x,
        AUC_time = AUC_time,
        auc_cal_method = auc_cal_method,
        risk_direction = "higher_is_worse"
      )
    })
  }


  list_input_data <- preprocess_previous_signature_survival_data(
    list_input_data,
    common_feature,
    context = "cal_auc_pre.prog.sig"
  )



  less.os.time <- max(list_input_data[[1]]$OS.time)

  for (i in names(list_input_data)) {
    print(i)
    if (max(list_input_data[[i]]$OS.time) < less.os.time) {
      less.os.time <- max(list_input_data[[i]]$OS.time)

      print(less.os.time)
    } else {
      less.os.time <- less.os.time
    }
  }

  print("Please wait a few minutes.")
  if (less.os.time > 365 * AUC_time) {
    model.name <- unique(sig.input$model)

    val_dd_list <- lapply(list_input_data, function(x) {
      x[, common_feature, drop = FALSE]
    })


    roc.table <- lapply(model.name, function(z) {
      coef.tab <- sig.input[sig.input$model == z, c("Coef", "symbol")]

      val_dd_list2 <- lapply(val_dd_list, function(x) {
        x[, c("OS.time", "OS", coef.tab$symbol), drop = FALSE]
      })


      rs <- lapply(val_dd_list2, function(x) {
        cbind(x[, 1:2],
          RS = apply(as.data.frame(x[, -c(1:2)]), 1, function(x) {
            x %*% coef.tab$Coef
          })
        )
      })

      roc.test <- returnRStoROC(rs.table.list = rs, AUC_time = AUC_time)




      return(roc.test)
    })
    names(roc.table) <- model.name


    return(roc.table)
  } else {
    stop(paste0(
      "AUC_time is outside the available follow-up window. The shortest ",
      "overall survival time in the provided cohorts is ", less.os.time,
      " days; choose AUC_time < ", round(less.os.time / 365, 3), " years."
    ), call. = FALSE)
  }
}
