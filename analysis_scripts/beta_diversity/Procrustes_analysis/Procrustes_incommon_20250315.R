library(phyloseq)
library(data.table)
library(ggplot2)
library(dplyr)
library(vegan)
library(ape)
library(readr)
library(tibble)
library(readxl)
library(stringr)
library(ggpubr)

#### ILO Cohort ####
# Load in metadata first:
# (Code from /restricted/projectnb/uh2-sebas/data/metagenomics/ILO_combined_cohort/processed_data/qc_analysis/Quality_Control_Metaphlan4_ILO_clean_09.30.2024.RMD)
meta.dir <-"/restricted/projectnb/uh2-sebas/data/metagenomics/ILO_combined_cohort/decoder/archived/"
phenofl <- "2024-11-18.cleaned.metag.list.csv"
#phenotypic data
pheno.data <- read_csv(paste0(meta.dir,phenofl) ) %>% 
  dplyr::mutate(Sample=stool_kit)
pheno.data$stool_kit <- as.character(pheno.data$stool_kit)
pheno.data$Sample <- as.character(pheno.data$Sample)


SAMPLE_phy <- sample_data(pheno.data %>% column_to_rownames(var="Sample"))

## Read in Bracken and MP4 phyloseq objects:

## Species

ILO.brkfile <- "/restricted/projectnb/uh2-sebas/data/metagenomics/ILO_combined_cohort/processed_data/data_library/bracken_renorm_ILO_phyloseq.R1.K3.10.09.2024.rds"
ILO.brk <- readRDS(ILO.brkfile)
sample_data(ILO.brk) <- SAMPLE_phy
ILO.brk.sp <- subset_samples(ILO.brk, !is.na(Age_at_FecalSample))
ILO.brk.sp.taxIDs <- rownames(otu_table(ILO.brk.sp))


ILO.mp4file <- "/restricted/projectnb/uh2-sebas/data/metagenomics/ILO_combined_cohort/processed_data/data_library/metaphlan4_renorm_ILO_phyloseq.09.30.2024.rds"
ILO.mp4 <- readRDS(ILO.mp4file)
sample_data(ILO.mp4) <- SAMPLE_phy
ILO.mp4.sp <- subset_samples(ILO.mp4, !is.na(Age_at_FecalSample))
ILO.mp4.sp.taxIDs <- rownames(otu_table(ILO.mp4.sp))

#Restrict to in common taxa and renormalize
otu_table(ILO.brk.sp) <- otu_table(ILO.brk.sp)[rownames(otu_table(ILO.brk.sp)) %in% intersect(ILO.brk.sp.taxIDs, ILO.mp4.sp.taxIDs),]
ILO.brk.sp <- transform_sample_counts(ILO.brk.sp, function(x) x / sum(x) )

otu_table(ILO.mp4.sp) <- otu_table(ILO.mp4.sp)[rownames(otu_table(ILO.mp4.sp)) %in% intersect(ILO.brk.sp.taxIDs, ILO.mp4.sp.taxIDs),]
ILO.mp4.sp <- transform_sample_counts(ILO.mp4.sp, function(x) x / sum(x) )

#Ordinations
pcoa_bracken_ILO <-  ordinate(ILO.brk.sp, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]

pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.sp, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)

test.result
## Genus

ILO.brk.ge <- tax_glom(ILO.brk.sp, "genus")
pcoa_bracken_ILO <-  ordinate(ILO.brk.ge, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]

ILO.mp4.ge <- tax_glom(ILO.mp4.sp, "Genus")
pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.ge, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)
test.result
## Family

ILO.brk.fa <- tax_glom(ILO.brk.sp, "family")
pcoa_bracken_ILO <-  ordinate(ILO.brk.fa, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]


ILO.mp4.fa <- tax_glom(ILO.mp4.sp, "Family")
pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.fa, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)
test.result

## Order

ILO.brk.or <- tax_glom(ILO.brk.sp, "order")
pcoa_bracken_ILO <-  ordinate(ILO.brk.or, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]


ILO.mp4.or <- tax_glom(ILO.mp4.sp, "Order")
pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.or, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)
test.result

## Class

ILO.brk.cl <- tax_glom(ILO.brk.sp, "class")
pcoa_bracken_ILO <-  ordinate(ILO.brk.cl, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]


ILO.mp4.cl <- tax_glom(ILO.mp4.sp, "Class")
pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.cl, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)

test.result


## Phylum

ILO.brk.ph <- tax_glom(ILO.brk.sp, "phylum")
pcoa_bracken_ILO <-  ordinate(ILO.brk.ph, "PCoA", "bray")
pcoa_bracken_ILO.vec <- as.data.frame(pcoa_bracken_ILO$vectors)
pcoa_bracken_ILO.vec <- pcoa_bracken_ILO.vec[order(rownames(pcoa_bracken_ILO.vec)),]


ILO.mp4.ph <- tax_glom(ILO.mp4.sp, "Phylum")
pcoa_metaphlan4_ILO <-  ordinate(ILO.mp4.ph, "PCoA", "bray")
pcoa_metaphlan4_ILO.vec <- as.data.frame(pcoa_metaphlan4_ILO$vectors)
pcoa_metaphlan4_ILO.vec <- pcoa_metaphlan4_ILO.vec[order(rownames(pcoa_metaphlan4_ILO.vec)),]

pro1 <- ade4::procuste(dfX = pcoa_bracken_ILO.vec,
                       dfY = pcoa_metaphlan4_ILO.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_ILO.vec,
                                    pcoa_metaphlan4_ILO.vec, nrepet = 999)

test.result

#### Xu et al Cohort ####
# Load in metadata first:
# (Code from restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/processed_data/qc_analysis/Quality_Control_Bracken_chinese_10.18.2024.RMD)
meta.dir <- "/restricted/projectnb/necs/Data_library/Metabolomics/resources/Xu_plasma_metabolite_analysis/"
sra.dir <- "/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/metadata/"
phenofl <- "phenotype.xlsx"

clinical.data <- read_excel(paste0(meta.dir, phenofl))

#add ethnicity variable
clinical.data$ethnicity <- "han_chinese"


#extract sample IDs and create family ID and group ID
sample.ID <- clinical.data %>% dplyr::select(1) %>% pull()

fam.group <- str_remove(sample.ID, "QCS-") %>%                       
  str_split("-") %>%
  unlist() 

fam.group <- matrix(fam.group, nrow=2)

#create family ID
fam.id <- fam.group[1,]

#create group ID
group.id <- fam.group[2,]
group.id[group.id == "1"] <- "Advanced Age"
group.id[group.id == "2" | group.id == "3"] <- "Offspring"


#add to clinical data
clinical.clean <- clinical.data %>%
  mutate(sample.ID = sample.ID, fam.id = fam.id, group = group.id) 
clinical.clean$`C-age (years)` <- as.integer(clinical.clean$`C-age (years)`)

#load sample Run names and sample SRA names
SRA.data <- read_csv(paste0(sra.dir,"SraRunInfo.csv"))

SRA.clean <- SRA.data %>% 
  dplyr::select(Run, SampleName) 

#sample ID and sample SRA names
sample.dict <- read_excel(paste0(sra.dir, "sample name-sequencing ID sheet.xlsx"))

#combine files
sample.info <- dplyr::inner_join(SRA.clean, sample.dict, by = c("SampleName"))

pheno.data <- dplyr::inner_join(sample.info, clinical.clean, by = c("sample.ID")) %>% 
  dplyr::select(-c(sample.ID)) %>%
  mutate(stool_kit = Run, Sample=Run) 
pheno.data$stool_kit <- as.character(pheno.data$stool_kit)
pheno.data$Sample <- as.character(pheno.data$Sample)

SAMPLE_phy <- sample_data(pheno.data %>% column_to_rownames(var="Sample"))

## Read in Bracken and MP4 phyloseq objects:

## Species plot

Xu.brkfile <- "/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/processed_data/data_library/bracken_renorm_chinese_phyloseq.R1.K3.10.22.2024.rds"
Xu.brk <- readRDS(Xu.brkfile)
sample_data(Xu.brk) <- SAMPLE_phy
Xu.brk.sp <- subset_samples(Xu.brk, !is.na(C.age..years.))
Xu.brk.sp.taxIDs <- rownames(otu_table(Xu.brk.sp))

Xu.mp4file <- "/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/processed_data/data_library/metaphlan4_renorm_chinese_phyloseq.10.23.2024.rds"
Xu.mp4 <- readRDS(Xu.mp4file)
sample_data(Xu.mp4) <- SAMPLE_phy
Xu.mp4.sp <- subset_samples(Xu.mp4, !is.na(C.age..years.))
Xu.mp4.sp.taxIDs <- rownames(otu_table(Xu.mp4.sp))

#Restrict to in common taxa and renormalize
otu_table(Xu.brk.sp) <- otu_table(Xu.brk.sp)[rownames(otu_table(Xu.brk.sp)) %in% intersect(Xu.brk.sp.taxIDs, Xu.mp4.sp.taxIDs),]
Xu.brk.sp <- transform_sample_counts(Xu.brk.sp, function(x) x / sum(x) )

otu_table(Xu.mp4.sp) <- otu_table(Xu.mp4.sp)[rownames(otu_table(Xu.mp4.sp)) %in% intersect(Xu.brk.sp.taxIDs, Xu.mp4.sp.taxIDs),]
Xu.mp4.sp <- transform_sample_counts(Xu.mp4.sp, function(x) x / sum(x) )

#Ordinate
pcoa_bracken_Xu <-  ordinate(Xu.brk.sp, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.sp, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                     pcoa_metaphlan4_Xu.vec, nrepet = 999)
test.result

## Genus

Xu.brk.ge <- tax_glom(Xu.brk.sp, "genus")
pcoa_bracken_Xu <-  ordinate(Xu.brk.ge, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)


Xu.mp4.ge <- tax_glom(Xu.mp4.sp, "Genus")
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.ge, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                                    pcoa_metaphlan4_Xu.vec, nrepet = 999)
test.result

## Family

Xu.brk.fa <- tax_glom(Xu.brk.sp, "family")
pcoa_bracken_Xu <-  ordinate(Xu.brk.fa, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)

Xu.mp4.fa <- tax_glom(Xu.mp4.sp, "Family")
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.fa, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                                    pcoa_metaphlan4_Xu.vec, nrepet = 999)

test.result


## Order

Xu.brk.or <- tax_glom(Xu.brk.sp, "order")
pcoa_bracken_Xu <-  ordinate(Xu.brk.or, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)


Xu.mp4.or <- tax_glom(Xu.mp4.sp, "Order")
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.or, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                                    pcoa_metaphlan4_Xu.vec, nrepet = 999)
test.result

## Class

Xu.brk.cl <- tax_glom(Xu.brk.sp, "class")
pcoa_bracken_Xu <-  ordinate(Xu.brk.cl, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)


Xu.mp4.cl <- tax_glom(Xu.mp4.sp, "Class")
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.cl, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                                    pcoa_metaphlan4_Xu.vec, nrepet = 999)
test.result

## Phylum

Xu.brk.ph <- tax_glom(Xu.brk.sp, "phylum")
pcoa_bracken_Xu <-  ordinate(Xu.brk.ph, "PCoA", "bray")
pcoa_bracken_Xu.vec <- as.data.frame(pcoa_bracken_Xu$vectors)

Xu.mp4.ph <- tax_glom(Xu.mp4.sp, "Phylum")
pcoa_metaphlan4_Xu <-  ordinate(Xu.mp4.ph, "PCoA", "bray")
pcoa_metaphlan4_Xu.vec <- as.data.frame(pcoa_metaphlan4_Xu$vectors)

pro1 <- ade4::procuste(dfX = pcoa_bracken_Xu.vec,
                       dfY = pcoa_metaphlan4_Xu.vec)
set.seed(42)
test.result <- ade4::procuste.rtest(pcoa_bracken_Xu.vec,
                                    pcoa_metaphlan4_Xu.vec, nrepet = 999)

test.result




