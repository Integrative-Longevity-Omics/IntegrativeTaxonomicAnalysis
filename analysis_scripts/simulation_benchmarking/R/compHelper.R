# Helper function to get stats from vector of p-values
evalFunc <- function(ps, alpha, iBoth, adjP = FALSE) {
  
  if(adjP) {
    ps <- p.adjust(ps, method = "BH")
  }
  
  sens <- ((ps < alpha) %*% iBoth) / sum(iBoth)
  spec <- ((ps > alpha) %*% !iBoth) / sum(!iBoth)
  prec <- ((ps < alpha) %*% iBoth) / sum((ps < alpha))
  if(sum(iBoth) > 0) {
    auc <- suppressMessages(pROC::auc(iBoth, ps))
    aucThre <- which.min(abs(attr(auc, "roc")$specificities - (1-alpha)))
    SensAlt <- attr(auc, "roc")$sensitivities[aucThre]
    pAlt <- sort(ps, decreasing = TRUE)[aucThre]
  } else {
    auc <- NaN
    SensAlt <- NaN
    pAlt <- NaN
  }
  return(c(sens, spec, prec, auc, SensAlt, pAlt))
}

#cResList <- lapply(c(0.05), alphaFunc, simI$metaPcorrMat, iBoth)


# Wrapper functions to get stats for given alpha
alphaFunc <- function(alpha, mPvals, iBoth) {
  
  # Nominal P-value
  evalMat1 <- t(apply(mPvals, 2, evalFunc, alpha, iBoth))
  evalMat1[is.nan(evalMat1)] <- 0
  colnames(evalMat1) <- c("sens", "spec", "prec", "auc", "sensAlt", "pAlt")
  evalMean1 <- colMeans(evalMat1)
  names(evalMean1) <- paste0(names(evalMean1), "_mean")
  res1 <- data.frame(
    t(evalMean1)
  )
  res1$p_type <- "nominal"
  
  # Format all results
  all1 <- as.data.frame(evalMat1)
  all1$p_type <- "nominal"
  all1$rep <- seq(nrow(all1))
  
  
  # FDR corrected
  evalMat2 <- t(apply(mPvals, 2, evalFunc, alpha, iBoth, TRUE))
  evalMat2[is.nan(evalMat2)] <- 0
  colnames(evalMat2) <- c("sens", "spec", "prec", "auc", "sensAlt", "pAlt")
  evalMean2 <- colMeans(evalMat2)
  names(evalMean2) <- paste0(names(evalMean2), "_mean")
  res2 <- data.frame(
    t(evalMean2)
  )
  res2$p_type <- "fdr"
  
  # Format all results
  all2 <- as.data.frame(evalMat2)
  all2$p_type <- "fdr"
  all2$rep <- seq(nrow(all2))
  
  # Compile results
  res <- rbind(res1, res2)
  res$alpha <- alpha
  
  all <- rbind(all1, all2)
  all$alpha <- alpha
  
  return(list(res = res,
              all = all))
  
}

#for naive approach
alphaFunc_naive <- function(alpha, mSignif_nominal, mSignif_fdr, iBoth) {
  
  #nominal p-values
  # Convert input to logical TRUE/FALSE
  mSignif_nominal <- mSignif_nominal == TRUE | mSignif_nominal == "yes"
  iBoth <- as.logical(iBoth)
  
  # Internal evaluation for one column
  evalOne <- function(pred, truth, alpha) {
    pred <- as.logical(pred)
    truth <- as.logical(truth)
    
    # Confusion matrix counts
    TP <- sum(pred & truth)
    FP <- sum(pred & !truth)
    TN <- sum(!pred & !truth)
    FN <- sum(!pred & truth)
    
    sens <- TP / (TP + FN)
    spec <- TN / (TN + FP)
    prec <- TP / (TP + FP)
    
    # AUC: convert binary predictions to numeric 0/1
    pred_num <- as.numeric(pred)
    
    if (sum(truth) > 0) {
      auc <- suppressMessages(pROC::auc(truth, pred_num))
      roc_obj <- attr(auc, "roc")
      aucThre <- which.min(abs(roc_obj$specificities - (1 - alpha)))
      SensAlt <- roc_obj$sensitivities[aucThre]
      pAlt <- sort(pred_num, decreasing = TRUE)[aucThre]
    } else {
      auc <- NaN
      SensAlt <- NaN
      pAlt <- NaN
    }
    
    return(c(sens, spec, prec, auc, SensAlt, pAlt))
  }
  
  # Apply to each column (simulation replicate)
  evalMat1 <- t(apply(mSignif_nominal, 2, evalOne, truth = iBoth, alpha = alpha))
  evalMat1[is.nan(evalMat1)] <- 0
  colnames(evalMat1) <- c("sens", "spec", "prec", "auc", "sensAlt", "pAlt")
  
  # Mean across columns
  evalMean1 <- colMeans(evalMat1)
  names(evalMean1) <- paste0(names(evalMean1), "_mean")
  res1 <- data.frame(t(evalMean1), p_type = "nominal", alpha = alpha)
  
  # Full per-column results
  all1 <- as.data.frame(evalMat1)
  all1$p_type <- "nominal"
  all1$rep <- seq(nrow(all1))
  
  
  #fdr p-values
  # Convert input to logical TRUE/FALSE
  mSignif_fdr <- mSignif_fdr == TRUE | mSignif_fdr == "yes"
  iBoth <- as.logical(iBoth)
  
  # Internal evaluation for one column
  evalTwo <- function(pred, truth, alpha) {
    pred <- as.logical(pred)
    truth <- as.logical(truth)
    
    # Confusion matrix counts
    TP <- sum(pred & truth)
    FP <- sum(pred & !truth)
    TN <- sum(!pred & !truth)
    FN <- sum(!pred & truth)
    
    sens <- TP / (TP + FN)
    spec <- TN / (TN + FP)
    prec <- TP / (TP + FP)
    
    # AUC: convert binary predictions to numeric 0/1
    pred_num <- as.numeric(pred)
    
    if (sum(truth) > 0) {
      auc <- suppressMessages(pROC::auc(truth, pred_num))
      roc_obj <- attr(auc, "roc")
      aucThre <- which.min(abs(roc_obj$specificities - (1 - alpha)))
      SensAlt <- roc_obj$sensitivities[aucThre]
      pAlt <- sort(pred_num, decreasing = TRUE)[aucThre]
    } else {
      auc <- NaN
      SensAlt <- NaN
      pAlt <- NaN
    }
    
    return(c(sens, spec, prec, auc, SensAlt, pAlt))
  }
  
  # Apply to each column (simulation replicate)
  evalMat2 <- t(apply(mSignif_fdr, 2, evalTwo, truth = iBoth, alpha = alpha))
  evalMat2[is.nan(evalMat2)] <- 0
  colnames(evalMat2) <- c("sens", "spec", "prec", "auc", "sensAlt", "pAlt")
  
  # Mean across columns
  evalMean2 <- colMeans(evalMat2)
  names(evalMean2) <- paste0(names(evalMean2), "_mean")
  res2 <- data.frame(t(evalMean2), p_type = "fdr", alpha = alpha)
  
  # Full per-column results
  all2 <- as.data.frame(evalMat2)
  all2$p_type <- "fdr"
  all2$rep <- seq(nrow(all2))
  
  # Compile results
  res <- rbind(res1, res2)
  res$alpha <- alpha
  
  all <- rbind(all1, all2)
  all$alpha <- alpha
  
  
  return(list(res = res, all = all))
}
