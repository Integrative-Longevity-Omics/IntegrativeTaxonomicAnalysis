# Function to run simulations based on a set of inputs parameters
simFunc <- function(corVal, alt = 0, propAlt = 0, both = 0, propBoth = 0, ncols = 2, nVal = 500) {
  
  # Create correlation matrix
  corMat <- matrix(corVal, ncols, ncols)
  diag(corMat) <- 1
  
  # Correct propAlt
  propAltCor <- propAlt * (1-propBoth)
  
  if(propAltCor > 0) {
    iAlt <- seq(nVal * propAltCor)
    nAlt <- sample(ncols-1, length(iAlt), replace = TRUE)
    wAlt <- lapply(nAlt, function(x) sample(ncols, x))
  } else {
    iAlt <- c()
  }
  
  if(propBoth > 0) {
    iBoth <- seq(nVal * propAltCor + nVal * propBoth)
    iBoth <- iBoth[!iBoth %in% iAlt]
  } else {
    iBoth <- c()
  }
  
  # Simulate test stastics and estimate correlation
  datSim <- mvrnorm(nVal, rep(0, ncol(corMat)), corMat)
  
  if(propAltCor > 0) {
    for(i in seq(length(wAlt))) {
      datSim[i, wAlt[[i]]] <- datSim[i, wAlt[[i]]] + alt
    }
  }
  
  # Add both
  datSim[iBoth,] <- datSim[iBoth,] + both
  
  corEst <- .tetCor(datSim)
  
  # Calculcate effective number of studies
  esEst <- .ENSest(corEst)
  
  # Convert test statistics to p-values
  datPvals <- pnorm(datSim, lower.tail = FALSE)
  
  # Estimate corrected meta-analysis p-values
  metaPcorr <- .adjMPest(datPvals, esEst)
  
  ## Province
  
  # Calculate Z variance
  varZ <- .VARest(corEst)
  
  # Calculate Z sums
  sumZ <- rowSums(datSim)
  
  # Calculate Province Z-score
  metaPprov <- .provBorEst(datSim, varZ)
  
  # Determine significance of shared features based on naive approach based on nominal p-value
  # return significance
  metaNaive_nominal <- compare_two_tests_nominal(datPvals,datSim,1,2,alpha=0.05)
  
  # Determine significance of shared features based on naive approach based on fdr
  # return significance
  metaNaive_fdr <- compare_two_tests_fdr(datPvals,datSim,1,2,alpha=0.05)
  
  return(list(
    metaPcorr = metaPcorr,
    metaPprov = metaPprov,
    metaNaive_nominal = metaNaive_nominal,
    metaNaive_fdr = metaNaive_fdr,
    corEst = mean(corEst[upper.tri(corEst)]),
    esEst = esEst,
    varEst = varZ
  ))
}

# Wrapper function for grid of parameters
simListFunction <- function(ncols, nfeats, Plist, cores = 1, sims) {
  
  mclapply(Plist, function(Pvec) {
    
    # Print progress
    wSim <- which(unlist(lapply(Plist, function(x) identical(x, Pvec))))
    print(paste(ncols, ":", sims, ":", wSim, "/", length(Plist)))
    
    # Get parameters
    cS <- Pvec[1]
    aS <- Pvec[2]
    pS <- Pvec[3]
    aB <- Pvec[4]
    pB <- Pvec[5]
    
    # Calculate True effective number of studies and variance
    tMat <- matrix(cS, ncols, ncols)
    diag(tMat) <- 1
    eValsTrue <- eigen(tMat)$values
    esTrue <- .ENSest(tMat)
    varTrue <- sum(tMat)
    
    # Run simulation
    PsimList <- do.call(rbind, mclapply(seq(sims), function(i) {
      simFunc(cS, aS, pS, aB, pB, ncols, nfeats)
    }, mc.cores = 1))
    
    # Compile results
    
    #Naive approach results
    
    ## Corrected meta-analysis p-value
    metaPcorrMat <- do.call(cbind, PsimList[, 1])
    metaPprovMat <- do.call(cbind, PsimList[, 2])
    metaNaive_nominal <- do.call(cbind, PsimList[, 3])
    metaNaive_fdr <- do.call(cbind, PsimList[, 4])
    corEstVec <- unlist(PsimList[, 5])
    esEstVec <-  unlist(PsimList[, 6])
    varEstVec <-  unlist(PsimList[, 7])
    
    # Return list of results
    return(list(
      metaPcorrMat = metaPcorrMat,
      metaPprovMat = metaPprovMat,
      metaNaive_nominal = metaNaive_nominal,
      metaNaive_fdr = metaNaive_fdr,
      corEstVec = corEstVec,
      esEstVec =  esEstVec,
      varEstVec = varEstVec,
      corTrue = cS,
      esTrue = esTrue,
      varTrue = varTrue,
      nTests = nfeats,
      altVal = aS,
      altProp = pS,
      bothVal = aB,
      bothProp = pB
    ))
  }, mc.cores = cores)
}