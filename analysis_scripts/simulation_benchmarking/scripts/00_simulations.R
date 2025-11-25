# Requirements

## Load project packages and directories
require(MASS)
require(parallel)
require(psych)
require(stats)

## Directories
baseDir <- ".."
resDir <- file.path(baseDir, "results")
rDir <- file.path(baseDir, "R")

## Set files to read and write

# Files to read in
  
# Files to write out
simResFile <- file.path(resDir, "simRes_NC_NF.rds")

## Read in helper functions
source(file.path(rDir, "simHelper.R"))
source(file.path(rDir, "runMethods.R"))

## Set parameters

# Parameters to test
Cvec <- seq(0, 1, by = 0.1)
Avec <- c(0, 2, 4)
PrVec <- c(0, 0.2)
Pgrid <- expand.grid(Cvec, Avec, PrVec, Avec, PrVec)
Pgrid <- Pgrid[!((rowSums(Pgrid[,c("Var2", "Var3")] != 0) == 1) | (rowSums(Pgrid[,c("Var4", "Var5")] != 0) == 1)),]
Pgrid <- Pgrid[order(-Pgrid$Var2),]
Plist <- as.list(as.data.frame(t(Pgrid)))
print(length(Plist))

# Number of features to simulation
Nvec <- c(300, 350, 400, 450)

# Set number of simulations to run
nSims <- 200

# Set number of cores to run
nCors <- 10

# Number of studies to simulate
nStud <- c(2)

## Run simulations for different parameters
invisible(
  lapply(nStud, function(nc) {
    lapply(Nvec, function(nf){
      RNGkind("L'Ecuyer-CMRG")
      set.seed(666)
      simRes <- simListFunction(nc, nf, Plist, nCors, nSims)
      names(simRes) <- NULL
      simResFileNC <- sub("NC", nc, simResFile)
      simResFileNC <- sub("NF", nf, simResFileNC)
      saveRDS(simRes, simResFileNC)
      gc()
    })
  })
)
