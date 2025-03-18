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
    #select(c(tax.name, contains("_bugs_list"))) %>%
    mutate(tax.id=tax.name) %>%
    select(-c(tax.name)) %>%
    column_to_rownames(var = "tax.id")
  dim(biom_prop)
  #useful::corner(biom_prop) #use this to visualize a corner of the data
  
  #rerun abundance table
  return(biom_prop)
  
}

reg_by_taxa <- function(indat, taxid, modformula, varlist){
  metadat <- as_tibble(indat) %>% 
    mutate(tax.ID = row.names(indat)) %>% 
    pivot_longer(!tax.ID) %>% 
    inner_join(pheno.data, by=c('name'='Sample.ID')) %>% 
    mutate(sample.id = name, value=as.numeric(value)) %>% 
    filter(tax.ID == taxid) %>%
    arrange(my.famID, sample.id) %>% 
    group_by(my.famID) %>%                       # Group by my.famID
    mutate(my.famID.idx = cur_group_id()) %>%      # Create a numeric index for each my.famID group
    ungroup() %>% 
    select(all_of(varlist), my.famID.idx) 
  # clean out missing
  metadat <- metadat[which(complete.cases(metadat)),]
  
  coeff.skew <- metadat %>% 
    select(c("tax.ID", "value"))  %>%
    group_by(tax.ID) %>%
    summarise_if(is.numeric, skewness) %>%
    rename("value"="skewness") 
  
  mod <- geese(modformula, id=my.famID.idx, corstr = "exchangeable", data=metadat)
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


p_relabudance_pheno <- function(indat=indat, tax.id=tax.id, phenovar){
  metadat.phylum <- as_tibble(indat) %>% 
    mutate(tax.ID = row.names(indat)) %>% 
    #pivot_longer(cols=contains("_bugs_list")) %>% 
    inner_join(pheno.data, by=c('name'='Sample')) %>% 
    mutate(sample.id = name, value=as.numeric(value)) %>% 
    filter(tax.ID == tax.id) %>%
    arrange(fam.id, sample.id)
  
  p <- ggplot(data=metadat.phylum, aes(x={{phenovar}}, y=value)) + geom_point()
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

pheno.data <- read.csv(opt$arg7)
# Sample.ID Age_at_FecalSample my.famID pn5_yrs_education
pheno.data$education <- ifelse(pheno.data$pn5_yrs_education %in% c(1,2,3,4), 1, 
                                 ifelse(pheno.data$pn5_yrs_education %in% c(7,8), 4, 
                                        ifelse(pheno.data$pn5_yrs_education == 5, 2,
                                               ifelse(pheno.data$pn5_yrs_education == 6, 3, pheno.data$pn5_yrs_education))))
pheno.data$Sample.ID <- as.character(pheno.data$stool_kit)

count.mat <- count_mat(taxadata = TAX, otudata=OTU, rank=opt$arg4)
prop.rel <- get_relabund(count.mat)
prop.rel.log <- log(prop.rel)


# Regression
taxid <- row.names(prop.rel.log)

f2 <- formula(value~Age_at_FecalSample+factor(registration_sex.factor))
res.f2 <- data.frame(t(sapply(taxid, reg_by_taxa, indat=prop.rel.log, modformula=f2, varlist=c("value","Age_at_FecalSample","registration_sex.factor","my.famID","tax.ID"))))

f4 <- formula(value~Age_at_FecalSample+factor(registration_sex.factor)+factor(education))
res.f4 <- data.frame(t(sapply(taxid, reg_by_taxa, indat=prop.rel.log, modformula=f4, varlist=c("value","Age_at_FecalSample","registration_sex.factor","my.famID","tax.ID", "education"))))

res.f2 <- out_parser(res.f2, opt$arg4)
res.f4 <- out_parser(res.f4, opt$arg4)


#export to csv
saveRDS(res.f2, 
        file=paste0(opt$arg2,paste0(opt$arg6, "_", opt$arg4, "_AgeGender_", opt$arg5, ".rds")))
saveRDS(res.f4, 
        file=paste0(opt$arg2,paste0(opt$arg6, "_", opt$arg4, "_AgeGenderEdu_", opt$arg5, ".rds")))

