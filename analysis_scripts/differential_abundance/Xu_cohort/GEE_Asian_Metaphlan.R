library(readxl)
library(dplyr)
library(ggplot2)
library(ape)
library(DESeq2)
library(stringr)
library(reshape2)
library(graphics)
library(optparse)

library(codebook)
library(future)
library(vtable)
library(useful)

library(vegan)
library(data.table)

# Data Cleaning
library(janitor)

# EDA
library(skimr)
library(DataExplorer)

# ggplot2 Helpers
library(ggplot2)
library(scales)
library(phyloseq)
library(taxonomizr)
library(tidyverse)
library(tidymodels)
library(egg)
library(gplots)
library(GGally)
library(ggcorrplot)
library(ggtext)
library(ggrepel)

#data transformations
library(sjmisc)

#PCA
library("FactoMineR")
library("factoextra")

#heatmap
library(pheatmap)

#MDMR
library(MDMR)

#Model
library(multilevelmod)
library(geepack)



# user defined functions
# calculate relative abundance/normalized abundance
get_relabund <- function(indat){
  relabund <- indat %>%
    as_tibble(rownames = "sample.ID") %>%
    janitor::adorn_percentages("col") %>%
    column_to_rownames(var = "sample.ID")
}

skewness <- function(x){
  (sum((x - mean(x))^3)/length(x))/
    ((sum((x - mean(x))^2)/length(x)))^(3/2)
}

kurtosis <- function(x){
  (sum((x - mean(x))^4)/length(x))/
    ((sum((x - mean(x))^2)/length(x)))^(4/2) - 3
}

# 
count_mat <- function(taxadata = biom_ilo_kraken.taxa, otudata=biom_ilo_kraken.otu, rank=rank){
  
  #convert taxonomic table to a tibble with rownames as a column called tax.ID
  #filter taxonomic table to remove taxa with no taxa classification in any level
  biom_taxdata <- as_tibble(taxadata, rownames = "tax.ID") %>% 
    unite("tax.name", kingdom:{{rank}}, remove = FALSE)
  
  #row number for taxa data after removing incomplete classification
  message(paste("taxa table rows:",nrow(biom_taxdata)))
  
  #convert otu table to a tibble with rownames as a column called tax.ID
  biom_otudata <- as_tibble(otudata, rownames = "tax.ID")
  #row number for otu data before removing incomplete classification
  message(paste("otu table rows:",nrow(biom_otudata)))
  
  #use an inner_join function from dplyr to join data from both table based on common tax ID
  biom_join <- inner_join(biom_taxdata, biom_otudata, by = "tax.ID")
  #row number for joined data
  message(paste("joined table rows:", nrow(biom_join))) #row number 11039 for the common taxa
  
  #remove rankings and transpose data for samples in rows and taxa IDs in columns
  #make sure the samples are set as rownames
  biom_prop <- biom_join %>% 
    group_by(tax.name) %>%
    summarise_if(is.numeric, sum) %>% 
    dplyr::select(c(tax.name, contains("SRR"))) %>%
    mutate(tax.id=tax.name) %>%
    dplyr::select(-c(tax.name)) %>%
    column_to_rownames(var = "tax.id")
  dim(biom_prop)
  #useful::corner(biom_prop) #use this to visualize a corner of the data
  
  #rerun abundance table
  return(biom_prop)
  
}

reg_by_taxa <- function(indat, taxid, modformula, varlist){
  metadat <- as_tibble(indat) %>% 
    mutate(tax.ID = row.names(indat)) %>% 
    pivot_longer(cols=contains("SRR")) %>% 
    inner_join(pheno.data, by=c('name'='Sample')) %>% 
    mutate(sample.id = name, value=as.numeric(value)) %>% 
    filter(tax.ID == taxid) %>%
    arrange(fam.id, sample.id) %>% 
    filter(!is.na(fam.id)) %>% 
    dplyr::select(all_of(varlist)) 
  # clean out missing
  metadat <- metadat[which(complete.cases(metadat)),]
  
  coeff.skew <- metadat %>% 
    dplyr::select(c("tax.ID", "value"))  %>%
    group_by(tax.ID) %>%
    summarise_if(is.numeric, skewness) %>%
    dplyr::rename(skewness=value) 
  
  mod <- geese(modformula, id=fam.id, corstr = "exchangeable", data=metadat)
  res <- summary(mod)$mean
  
  res <- res %>% 
    rownames_to_column(var="term") %>% 
    filter(term != "(Intercept)") %>% 
    dplyr::summarize(across(everything(), list)) %>% 
    filter(term != ",") %>%
    cbind(coeff.skew[,c("skewness")])
  #print(res)
  return(res)
}


manhplot <- function(est_dat, taxadata, rank, p, sig){
  biom_taxdata <- as_tibble(taxadata, rownames = "tax.ID") %>% 
    unite("tax.name", kingdom:species, remove = FALSE) %>% 
    unite("rank.name", kingdom:{{rank}}, remove=F)
  
  est_dat <- est_dat %>% 
    rownames_to_column(var="tax.name") %>% 
    left_join(biom_taxdata, by = "tax.name") %>%
    group_by(rank.name) %>% 
    add_count() %>% 
    arrange(n, rank.name, tax.name)
  
  est_dat <- est_dat  %>% 
    rowid_to_column("idx_cum")
  
  axis_set <- est_dat %>% 
    group_by(rank.name) %>% 
    summarize(center = mean(idx_cum)) %>% 
    add_count() %>%
    arrange(n, rank.name)
  axis_set$rank.name <- factor(axis_set$rank.name, levels=axis_set$rank.name)
  axis_set$rank.name <- stringr::str_wrap(axis_set$rank.name, 10)
  
  ylim <- ceiling(-log10(min(est_dat$Pval))) + 2
  
  p <- ggplot(est_dat, aes(x = idx_cum, y = -log10({{p}}), 
                           color = as_factor(rank.name), size = -log10({{p}}))) +
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.75) +
    facet_grid(cols=vars(rank.name),
               space="free_x",
               scales="free_x",
               switch="x") + 
    geom_text_repel(aes(label=ifelse(-log10({{p}}) >= -log10(sig), stringr::str_wrap(as.character(species.x), 10), '')), box.padding=0.5, segment.color='grey50', size=2, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) + 
    scale_x_continuous(label = axis_set$rank.name, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), unique(length(axis_set$rank.name)))) +
    scale_size_continuous(range = c(0.5,1)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5),
      panel.spacing = unit(.2, "cm"),
      strip.text = element_blank()
    )
  return(p)
}

circular_plot <- function(est_dat, taxadata, rank, pval, sig){
  biom_taxdata <- as_tibble(taxadata, rownames = "tax.ID") %>% 
    unite("tax.name", kingdom:species, remove = FALSE) %>% 
    unite("rank.name", kingdom:{{rank}}, remove=F)
  
  est_dat <- est_dat %>% 
    rownames_to_column(var="tax.name") %>% 
    left_join(biom_taxdata, by = "tax.name") %>%
    arrange(desc(tax.name)) %>% 
    mutate(logp = -log10({{pval}}))
  
  
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(factor(est_dat$rank.name)), ncol(est_dat)) )
  colnames(to_add) <- colnames(est_dat)
  to_add$rank.name <- rep(levels(factor(est_dat$rank.name)), each=empty_bar)
  est_dat <- rbind(est_dat, to_add)
  est_dat <- est_dat %>% 
    arrange(rank.name) %>% 
    mutate(id = row_number())
  
  
  # Get the name and the y position of each label
  label_data <- est_dat
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data <- est_dat %>% 
    group_by(rank.name) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end))) 
  
  base_data$id <- seq(1, nrow(base_data))
  
  number_of_bar <- nrow(base_data)
  angle <- 90 - 360 * (base_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  base_data$hjust <- ifelse(angle < -90, 1, 0)
  base_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  # Make the plot
  p <- ggplot(est_dat, aes(x=as.factor(id), y=logp*3, fill=rank.name)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=logp*3, fill=rank.name), stat="identity", alpha=0.5) +
    
    # Add a val=6/4/2 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = -log10(sig)*3, xend = start, yend = -log10(sig)*3), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    
    geom_bar(aes(x=as.factor(id), y=logp*3, fill=rank.name), stat="identity", alpha=0.5) +
    ylim(-6,30) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data=base_data, aes(x=title, y=15, label=rank.name, hjust=hjust, colour=rank.name), fontface="bold",alpha=1, size=2, angle= base_data$angle, inherit.aes = FALSE ) 
  
  p
}

out_parser <- function(indat, taxa){
  n.test <- nrow(indat)
  #sig <- 0.05/n.test
  sig <- 0.25
  indat <- data.frame(indat)
  
  s_vals <- c(".*k__",".*p__",".*c__",".*o__",".*f__",".*g__",".*s__")
  names(s_vals) <- c("kingdom","phylum","class","order","family","genus","species")
  taxapattern <- s_vals[taxa]
  
  indat$Est <- as.numeric(lapply(lapply(indat$estimate, '[[', 1), '[[',1))
  indat$Pval <- as.numeric(lapply(lapply(indat$p, '[[', 1), '[[',1))
  indat$Effect_Estimate <- factor(ifelse(indat$Est >= 0, "Increasing", "Decreasing"))
  
  indat$p.adj <- p.adjust(indat$Pval, method = "BH")
  indat$p.adj <- ifelse(indat$p.adj == 0, indat$p.adj + 0.000001, indat$p.adj)
  indat$Skew <- as.numeric(lapply(indat$skewness, '[[', 1))
  
  indat$taxaname <- sub(taxapattern, "", row.names(indat))
  
  return(indat[,c("Est", "Pval", "Effect_Estimate", "p.adj", "taxaname")])
}

p.adj.mod <- function (p, method = p.adjust.methods, n = length(p)) 
{
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  else p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp)
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  p0
}

p_relabudance_pheno <- function(indat=indat, tax.id=tax.id, phenovar){
  metadat.phylum <- as_tibble(indat) %>% 
    mutate(tax.ID = row.names(indat)) %>% 
    pivot_longer(cols=contains("SRR")) %>% 
    inner_join(pheno.data, by=c('name'='Sample')) %>% 
    mutate(sample.id = name, value=as.numeric(value)) %>% 
    filter(tax.ID == tax.id) %>%
    arrange(fam.id, sample.id)
  
  p <- ggplot(data=metadat.phylum, aes(x={{phenovar}}, y=value)) + geom_point()
  p
}

# Define the arguments
option_list <- list(
  make_option(c("-f", "--arg1"), type="character", default=NULL, help="the path to working files", metavar="character"),
  make_option(c("-o", "--arg2"), type="character", default=NULL, help="the directory of output files", metavar="character"),
  make_option(c("-b", "--arg3"), type="character", default=NULL, help="the directory of biom files", metavar="character"),
  make_option(c("-t", "--arg4"), type="character", default=NULL, help="the level of taxa", metavar="character"),
  make_option(c("-n", "--arg5"), type="character", default=NULL, help="the name of output files", metavar="character"),
  make_option(c("-m", "--arg6"), type="character", default=NULL, help="Regression method (Linear/GEE)", metavar="character"),
  make_option(c("-p", "--arg7"), type="character", default=NULL, help="the path to pheno files", metavar="character")
)

# Parse the arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



# Read in data
filtered_biom <- readRDS(paste0(opt$arg1,opt$arg3))

OTU <- filtered_biom@otu_table@.Data
TAX <- filtered_biom@tax_table@.Data 
colnames(TAX) <- tolower(colnames(TAX))

# 1. Read in the first file with the sheet "Table S1 Metadata"
metadata <- read_excel(opt$arg7,
                       sheet = "Table S1 Metadata", skip=1)
namelist <- make.names(names(metadata))
namelist[1] <- "Subject.No"
names(metadata) <- namelist
print(names(metadata))

# Rename the first column


# 2. Read in the second file with the sample name and sequencing ID sheet
sample_mapping <- read_excel("/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/metadata/sample name-sequencing ID sheet.xlsx", sheet="rename")

# 3. Read in the SraRunInfo.csv file with LibraryName and Run columns
sra_info <- read.csv("/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/metadata/SraRunInfo.csv")

# Perform the mappings based on the instructions
# Match Subject.No with sample.ID
merged_data <- metadata %>%
  dplyr::rename(sample.ID = Subject.No) %>%
  left_join(sample_mapping, by = "sample.ID")

# Match SampleName with LibraryName
pheno.data <- merged_data %>%
  left_join(sra_info %>% select(Run, LibraryName), by = c("SampleName" = "LibraryName")) %>% 
  select(-c(SampleName)) %>%
  dplyr::rename(Subject.idx = sample.ID, Sample=Run)

rm(merged_data)

pheno.data <- pheno.data %>%
  # Extract the middle three digits using regular expressions
  mutate(family_code = str_extract(Subject.idx, "(?<=QCS-)[0-9]{3}(?=-)")) %>%
  # Create a numeric index (fam.id) based on the order of unique middle three digits
  mutate(fam.id = as.numeric(factor(family_code))) %>%
  # Drop the intermediate 'family_code' column if no longer needed
  select(-family_code)

pheno.data$C.age..years. <- as.numeric(pheno.data$C.age..years.)
pheno.data$ethnicity <- "han_chinese"

count.mat <- count_mat(taxadata = TAX, otudata=OTU, rank=opt$arg4)

prop.rel <- get_relabund(count.mat)
prop.rel.log <- log(prop.rel)


# Regression
taxid <- row.names(prop.rel.log)


f2 <- formula(value~C.age..years.+factor(C.gender.1.male.2.female.))
res.f2 <- data.frame(t(sapply(taxid, reg_by_taxa, indat=prop.rel.log, modformula=f2, varlist=c("value","C.age..years.","C.gender.1.male.2.female.","fam.id","tax.ID"))))

f4 <- formula(value~C.age..years.+factor(C.gender.1.male.2.female.)+as.factor(E.years.of.education))
res.f4 <- data.frame(t(sapply(taxid, reg_by_taxa, indat=prop.rel.log, modformula=f4, varlist=c("value","C.age..years.","C.gender.1.male.2.female.",
                                                                                               "E.years.of.education","fam.id","tax.ID"))))
res.f2 <- out_parser(res.f2, opt$arg4)
res.f4 <- out_parser(res.f4, opt$arg4)

#export to csv
saveRDS(res.f2, 
        file=paste0(opt$arg2,paste0(opt$arg6, "_", opt$arg4, "_AgeGender_", opt$arg5, ".rds")))
saveRDS(res.f4, 
        file=paste0(opt$arg2,paste0(opt$arg6, "_", opt$arg4, "_AgeGenderEdu_", opt$arg5, ".rds")))
