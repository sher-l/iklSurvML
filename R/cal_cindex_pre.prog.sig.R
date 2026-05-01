#' Calculating the C index of the signatures from the previous paper
#'
#' @param use_your_own_collected_sig Whether to use your own collected signatures. T or F.
#' @param collected_sig_table If use_your_own_collected_sig is set as T, you should provide a data frame containing the information of the signatures. The column names of the data frame are "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol". For example, 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'. The 'model' consists of the first name of the first author and the PMID. The 'PMID' is from the paper. The 'Cancer' uses abbreviations in the format of the TCGA. 'Author' is the first name of the first author of the paper. 'Coef' is the coefficient of the variable in the signature. 'symbol' is the variable in the signature. If use_your_own_collected_sig is F, you don't need to provide this data.
#' @param type.sig  If the use_your_own_collected_sig is F, here we collected some signatures of the LGG, GBM, and Glioma. You can choose c('Glioma', 'LGG', 'GBM'), c('Glioma'),c('GBM'), c('Glioma', 'LGG'), and some other combination for the signatures you want.
#' @param list_input_data A list of the cohorts. Column names are 'ID', 'OS.time', 'OS', and the other variables. OS.time (Day). OS(1/0).
#'
#' @return A list of the C index score of the signature for each cohort you provide.
#' @export
#'
cal_cindex_pre.prog.sig <- function(use_your_own_collected_sig, # 是否使用您自己收集的signature， T or F
                                    collected_sig_table, # 列名分别为
                                    # "model"  "PMID"   "Cancer" "Author" "Coef"   "symbol"
                                    # 'Chen.33591634' '33591634''LGG' 'Chen' '0.7426' 'FGF7'
                                    type.sig, ### prognostic signature 的类型，c('Glioma','LGG','GBM')， c('Glioma'),c('Glioma','LGG')
                                    list_input_data # list of the cohorts(ID,OS.time, OS····)经过了log2（x+1）转化
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

  list_input_data <- preprocess_previous_signature_survival_data(
    list_input_data,
    common_feature,
    context = "cal_cindex_pre.prog.sig"
  )



  model.name <- unique(sig.input$model)

  val_dd_list <- lapply(list_input_data, function(x) {
    x[, common_feature, drop = FALSE]
  })


  cc.table <- lapply(model.name, function(z) {
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

    cc <- data.frame(Cindex = sapply(rs, function(x) {
      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
    })) %>%
      rownames_to_column("ID")



    return(cc)
  })
  names(cc.table) <- model.name


  return(cc.table)
}
