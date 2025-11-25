# Function to run AdjMaxP
runAdjMaxP <- function(
    pval_mat # Rows:Shared Features, Cols:Individual Studies
) {
  
  # Convert p-values to z-scores
  z_mat <- qnorm(pval_mat, lower.tail = FALSE)
  
  # Get tetrachoric correlation estimate
  corEst <- .tetCor(z_mat)
  
  # Calculate ENS
  ensEst <- .ENSest(corEst)
  
  # Calculate adjMaxPs
  metaPcorr <- .adjMPest(pval_mat, ensEst)
  
  # Return list of results
  outList <- list(
    adjMaxP_pvals = metaPcorr, # The pvalue estimates from AdjMax
    eff_no_studies = ensEst, # The estimate of the effective number of studies
    cor_est = corEst, # Inter-study correlation matrix
    z_mat = z_mat # Probit transformed p-values, z-score matrix (for plotting)
  )
  
}

#Function to run naive intersection approach
compare_two_tests_nominal <- function(pmat, tmat, test1, test2, alpha = 0.05) {
  # Extract columns
  p1 <- pmat[, test1]
  p2 <- pmat[, test2]
  t1 <- tmat[, test1]
  t2 <- tmat[, test2]
  
  # Significance indicators
  sig1 <- p1 < alpha
  sig2 <- p2 < alpha
  
  # Same direction (both positive or both negative)
  same_dir <- sign(t1) == sign(t2)
  
  # Combined condition
  shared_sig_same_dir <- sig1 & sig2 & same_dir
  
  # Return data frame
  # data.frame(
  #   feature = rownames(pmat),
  #   p_test1 = p1,
  #   p_test2 = p2,
  #   t_test1 = t1,
  #   t_test2 = t2,
  #   sig_test1 = sig1,
  #   sig_test2 = sig2,
  #   same_direction = same_dir,
  #   shared_sig_same_dir = ifelse(shared_sig_same_dir, "yes", "no"),
  #   row.names = NULL
  # )
  return(shared_sig_same_dir)
}

#Function to run naive intersection approach
compare_two_tests_fdr <- function(pmat, tmat, test1, test2, alpha = 0.05) {
  # Extract columns
  p1 <- pmat[, test1]
  p2 <- pmat[, test2]
  t1 <- tmat[, test1]
  t2 <- tmat[, test2]
  
  #fdr adjusted
  p1 <- p.adjust(p1, method = "BH")
  p2 <- p.adjust(p2, method = "BH")
  
  # Significance indicators
  sig1 <- p1 < alpha
  sig2 <- p2 < alpha
  
  # Same direction (both positive or both negative)
  same_dir <- sign(t1) == sign(t2)
  
  # Combined condition
  shared_sig_same_dir <- sig1 & sig2 & same_dir
  
  # Return data frame
  # data.frame(
  #   feature = rownames(pmat),
  #   p_test1 = p1,
  #   p_test2 = p2,
  #   t_test1 = t1,
  #   t_test2 = t2,
  #   sig_test1 = sig1,
  #   sig_test2 = sig2,
  #   same_direction = same_dir,
  #   shared_sig_same_dir = ifelse(shared_sig_same_dir, "yes", "no"),
  #   row.names = NULL
  # )
  return(shared_sig_same_dir)
}

# Function to run AdjMaxP
runProvBor <- function(
    pval_mat # Rows:Shared Features, Cols:Individual Studies
) {
  
  # Convert p-values to z-scores
  z_mat <- qnorm(pval_mat, lower.tail = FALSE)
  
  # Get tetrachoric correlation estimate
  corEst <- .tetCor(z_mat)
  
  # Return list of results
  outList <- list(
    provbor_pvals = metaPprov, # The pvalue estimates from province borecks
    mod_variant = varZ, # The estimate of the effective number of studies
    cor_est = corEst, # Inter-study correlation matrix
    z_mat = z_mat # Probit transformed p-values, z-score matrix (for plotting)
  )
  
}


.tetCor <- function(mat) {
  # Estimate tetrachoric correlation
  corE <- suppressMessages(tetrachoric(mat > 0, correct=FALSE, smooth = TRUE))$rho
  corE[corE < 0] <- 0
  return(corE)
}

# Calculate effective number of tests
.ENSest <- function(corE) {
  eVals <- eigen(corE)$values-1
  eVals[eVals < 0] <- 0
  ensE <- length(eVals) - sum(eVals)
  return(ensE)
}

# Calculate adjMaxPs
.adjMPest <- function(pmat, ensE) {
  
    # Get maximum P values
    maxPs <- apply(pmat, 1, max)
    
    # Calculate corrected p-values
    metaPCs <- pbeta(maxPs, ensE, 1) 
    
    return(metaPCs)
    
}

# Calculate Z variance
.VARest <- function(corE) {
  
  varZ <- sum(corE)
  
  return(varZ)
  
}

# Calculate Prov-Borecki p-value
.provBorEst <- function(zmat, varz) {
  
  # Calculate Z sums
  sumZ <- rowSums(zmat)
  
  # Calculate Province Z-score
  metaPprov <- pnorm(sumZ, 0, sqrt(varz), lower.tail = FALSE)
  
  return(metaPprov)
  
}

