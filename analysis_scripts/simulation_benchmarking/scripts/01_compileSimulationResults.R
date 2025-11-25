# Requirements

## Load project packages and directories

## Load packages and directories
require(MASS)
require(parallel)
require(psych)
require(stats)
require(pROC)
require(dplyr)
require(tidyr)
require(openxlsx)

## Directories
baseDir <- ".."
rDir <- file.path(baseDir, "R")
resDir <- file.path(baseDir, "results")
compDir <- file.path(resDir, "compiled")
tabDir <- file.path(resDir, "tables")

## Parameters
nCors <- 12

## Set files to read and write

# Files to read in
  
## Simulation results
simResFileList <- list.files(resDir, pattern = "simRes_")

# Files to write out
simCompFile <- file.path(compDir, "compiled_performance.rds")
simSumFile <- file.path(tabDir, "SimulationResults.xlsx")

## Read in helper functions
source(file.path(rDir, "compHelper.R"))

### Compile results
compiledPerf <- do.call(rbind, lapply(simResFileList, function(simResFile) {
  
  print(simResFile)
  
  simRes <- readRDS(file.path(resDir, simResFile))
  Nstudies <- unlist(strsplit(simResFile, "_"))[2]
  
  ## Get stats for each model
  simResAggList <- do.call(rbind, mclapply(simRes, function(simI) {
    
    # Pvalues under the alternative hypothesis
    altPropCor <- simI$altProp*(1-simI$bothProp)
    
    if(altPropCor > 0) {
      iAlt <- seq(simI$nTests)[seq(simI$nTests * altPropCor)]
    } else {
      iAlt <- c()
    }
    
    if(simI$bothProp > 0) {
      iBoth <- seq(simI$nTests) %in% seq(simI$nTests * altPropCor + simI$nTests * simI$bothProp)
      iBoth[iAlt] <- FALSE
    } else {
      iBoth <- rep(FALSE, simI$nTests)
    }
    
    # Corrected model
    cResList <- lapply(0.05, alphaFunc, simI$metaPcorrMat, iBoth)
    
    cRes <- do.call(rbind, lapply(cResList, function(x) x$res))
    cRes$meth <- "amp"
    
    cResAll <- do.call(rbind, lapply(cResList, function(x) x$all))
    cResAll$meth <- "amp"
    
    # Province model
    pResList <- lapply(0.05, alphaFunc, simI$metaPprovMat, iBoth)
    
    pRes <- do.call(rbind, lapply(pResList, function(x) x$res))
    pRes$meth <- "pb"
    
    pResAll <- do.call(rbind, lapply(pResList, function(x) x$all))
    pResAll$meth <- "pb"
    
    # Naive approach
    nResList <- lapply(0.05, alphaFunc_naive, simI$metaNaive_nominal, simI$metaNaive_fdr,iBoth)
    
    nRes <- do.call(rbind, lapply(nResList, function(x) x$res))
    nRes$meth <- "naive"
    
    nResAll <- do.call(rbind, lapply(nResList, function(x) x$all))
    nResAll$meth <- "naive"
    
    # Combine results
    eRes <- rbind(cRes, pRes, nRes)
    eRes$p_type <- as.factor(eRes$p_type)
    eRes$alpha <- as.factor(eRes$alpha)
    eRes$meth <- as.factor(eRes$meth)
    
    eResAll <- rbind(cResAll, pResAll, nResAll)
    eResAll$p_type <- as.factor(eResAll$p_type)
    eResAll$alpha <- as.factor(eResAll$alpha)
    eResAll$meth <- as.factor(eResAll$meth)
    
    ## Return data frame of results
    simIdf <- data.frame(
      nStudies = Nstudies,
      nTest = simI$nTest,
      alt = simI$altVal,
      altProp = simI$altProp,
      both = simI$bothVal,
      bothProp = simI$bothProp,
      corP = simI$corTrue,
      et = simI$esTrue,
      varP = simI$varTrue,
      corEmean = mean(simI$corEstVec),
      etEmean = mean(simI$esEstVec),
      varEmean = mean(simI$varEstVec),
      eRes
    )
    
    simIdfALL <- data.frame(
      nStudies = Nstudies,
      nTest = simI$nTest,
      alt = simI$altVal,
      altProp = simI$altProp,
      both = simI$bothVal,
      bothProp = simI$bothProp,
      corP = simI$corTrue,
      et = simI$esTrue,
      varP = simI$varTrue,
      corEmean = simI$corEstVec,
      etEmean = simI$esEstVec,
      varEmean = simI$varEstVec,
      eResAll
    )
    
    return(list(simIdf,
                simIdfALL))
    
  }, mc.cores = nCors))
  
  saveRDS(simResAggList, paste0(resDir, "/simResAggList.rds"))
  
  rm(simRes)
  gc()
  
  simResAggSub <- do.call(rbind, simResAggList[,1])
  rownames(simResAggSub) <- NULL
  
  simResAggAllSub <- do.call(rbind, simResAggList[,2])
  rownames(simResAggAllSub) <- NULL
  
  return(list(
    simResAggSub,
    simResAggAllSub
  ))
  
}))

saveRDS(compiledPerf,paste0(resDir, "/simResAggAllSel.rds"))

# Aggregate all results
simResAgg <- do.call(rbind, compiledPerf[,1])
rownames(simResAgg) <- NULL

simResAggAll <- do.call(rbind, compiledPerf[,2])
rownames(simResAggAll) <- NULL

# Combine parameter values

simResGroups <- unique(simResAgg[, c("bothProp", "both", "altProp", "alt")])
simResGroups$signalCons <- paste0("Conserved Features\nProportion=", simResGroups$bothProp,
                              ", Signal=", simResGroups$both)
simResGroups$signalCons[simResGroups$signalCons == "Conserved Features\nProportion=0, Signal=0"] <- "No Conserved Features"
consLevels <- c("No Conserved Features", "Conserved Features\nProportion=0.2, Signal=2")
simResGroups <- simResGroups[simResGroups$signalCons %in% consLevels,]
simResGroups$signalCons <- factor(simResGroups$signalCons, levels = consLevels)

simResGroups$signalUnCons <- paste0("Unconserved Features\nProportion=", simResGroups$altProp,
                              ", Signal=", simResGroups$alt)
simResGroups$signalUnCons[simResGroups$signalUnCons == "Unconserved Features\nProportion=0, Signal=0"] <- "No Unconserved Features"
UnConsLevels <- c("No Unconserved Features", "Unconserved Features\nProportion=0.2, Signal=2", "Unconserved Features\nProportion=0.2, Signal=4")
simResGroups <- simResGroups[simResGroups$signalUnCons %in% UnConsLevels,]
simResGroups$signalUnCons <- factor(simResGroups$signalUnCons, levels = UnConsLevels)

simResAggSel <- merge(simResGroups, simResAgg)
simResAggAllSel <- merge(simResGroups, simResAggAll)

# Format data for plotting
simResAggAllSel$corPfac <- as.factor(simResAggAllSel$corP)
simResAggAllSel$corPfac <- recode(simResAggAllSel$corPfac, "0" = "0.0", "1" = "1.0")

# Save results
saveRDS(simResAggAllSel, simCompFile)

# Format and save simulation results
simResAggOut <- simResAggSel %>%
  filter(alpha == 0.05) %>%
  dplyr::select(-c(signalCons, signalUnCons, etEmean, varEmean, et, varP, alpha)) %>%
  mutate(p_type = as.character(p_type))
colnames(simResAggOut) <- sub("_mean|mean", "", colnames(simResAggOut))

simResAggOut <- simResAggOut %>%
  gather("measure", "mean_val", c(corE, spec, sens, prec, auc, sensAlt, pAlt))

simResAggOut <- simResAggOut %>%
  pivot_wider(id_cols = c(nStudies, nTest, alt, both, altProp, bothProp, p_type, measure),
              names_from = c(meth, corP),
              names_sort = TRUE,
              values_from = mean_val) %>%
  filter(!((p_type == "fdr" & measure == "corE") |
             (p_type == "fdr" & measure == "auc") |
             (p_type == "fdr" & measure == "sensAlt") |
             (p_type == "fdr" & measure == "pAlt") |
             (both == 0 & measure %in% c("auc", "prec", "sens", "sensAlt", "pAlt")) |
             (altProp == 0.3) |
             (bothProp == 0.3) |
             (both == 4)))

simResAggOut$p_type[simResAggOut$measure == "corE"] <- "corE"
simResAggOut$p_type[simResAggOut$measure == "auc"] <- "auc"
simResAggOut$p_type[simResAggOut$measure == "sensAlt"] <- "Sens. (Spec. = 0.95)"

simResAggOut$measure <- recode_factor(simResAggOut$measure,
                                      "corE" = "Cor. Est.",
                                      "auc" = "AUC",
                                      "spec" = "Specificity",
                                      "sens" = "Sensitivity",
                                      "sensAlt" = "Sens. (Spec. = 0.95)",
                                      "pAlt" = "P-value (Spec. = 0.95)",
                                      "prec" = "Precision")

simResAggOut$p_type <- recode_factor(simResAggOut$p_type,
                                     "corE" = "None",
                                     "auc" = "None",
                                     "nominal" = "P-value < 0.05",
                                     "fdr" = "FDR Q-value < 0.05")

simResAggOut <- simResAggOut  %>%
  arrange(nStudies, nTest, alt, altProp, both, bothProp, p_type, measure)

# Write results to excel book
invisible(suppressWarnings(file.remove(simSumFile)))
wb <- createWorkbook()

invisible(lapply(unique(simResAggOut$nStudies), function(n) {

  sheetName <- paste0("Studies", n)
  simResAggOutN <- simResAggOut %>%
    filter(nStudies == n)

  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, simResAggOutN, rowNames =FALSE)

}))

saveWorkbook(wb, simSumFile, overwrite = TRUE)

