require(MASS)
require(psych)

# Function to run AdjMaxP
runAdjMaxP <- function(
  pval_mat # Rows:Shared Features, Cols:Individual Studies
) {
  
  # Convert p-values to z-scores
  z_mat <- qnorm(pval_mat, lower.tail = FALSE)
  
  # Estimate tetrachoric correlation
  corEst <- suppressMessages(tetrachoric(z_mat > 0, correct=FALSE, smooth = TRUE,))$rho
  corEst[corEst < 0] <- 0
  
  # Calculate effective number of tests
  eVals <- eigen(corEst)$values-1
  eVals[eVals < 0] <- 0
  ensEst <- length(eVals) - sum(eVals)
  
  # Calculate maximum p-values
  maxPs <- apply(pval_mat, 1, max)
  
  # Calculate corrected p-values
  metaPcorr <- pbeta(maxPs, ensEst, 1)
  
  # Return list of results
  outList <- list(
    adjMaxP_pvals = metaPcorr, # The pvalue estimates from AdjMax
    eff_no_studies = ensEst, # The estimate of the effective number of studies
    cor_est = corEst, # Inter-study correlation matrix
    z_mat = z_mat # Probit transformed p-values, z-score matrix (for plotting)
  )
  
}

# Example

## Simulate p-values from two correlated studies

### Toy 2x2 correlation matrix
corSim <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

### Simulated p-values from statistical tests
zSim <- mvrnorm(500, rep(0, ncol(corSim)), corSim)
pSim <- pnorm(zSim, lower.tail = FALSE)

dim(pSim)
# [1] 200   2

### Run AdjMaxP
res <- runAdjMaxP(pSim)

### AdjMaxP p-values
hist(res$adjMaxP_pvals)
#### This null example will be approximately uniform between 0 and 1

### Effective number of studies
res$eff_no_studies
#### Will be ~1.5 given corSim and 2 studies

### Tetrachoric correlation matrix
res$cor_est
#### Will be similar to corSim

### Probit transformed p-values to z-scores
plot(res$z_mat)
#### Will look correlated

