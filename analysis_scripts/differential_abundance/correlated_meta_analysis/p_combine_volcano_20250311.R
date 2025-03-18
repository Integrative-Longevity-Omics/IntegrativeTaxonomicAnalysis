library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tidyr)
library(dplyr)

setwd("/restricted/projectnb/uh2-sebas/analysis/metagenomics/meg_analyses")
source("adjMaxP_codePlusExample.R")
indir <- "/restricted/projectnb/uh2-sebas/analysis/metagenomics/ye_analyses/Output/"
outdir <- "volcano_plots"
#dir.create(outdir)

### ILO Cohort ###

reg_results_ILO <- readxl::read_excel(file.path(indir, "Common_taxa_regression.xlsx"), sheet = "R1K3")
colnames(reg_results_ILO)

pvals <- as.matrix(reg_results_ILO %>% dplyr::select(Pval_Bracken, Pval_Metaphlan))
dim(pvals)

### Run AdjMaxP
res <- runAdjMaxP(pvals)

### AdjMaxP p-values
hist(res$adjMaxP_pvals, breaks=100)
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

# Investigate odd pattern
hist(reg_results_ILO$Pval_Bracken, breaks = 1000)
hist(reg_results_ILO$Pval_Metaphlan, breaks = 1000)

fdrq <- p.adjust(res$adjMaxP_pvals, method="BH")

reg_results_ILO$AdjMaxP <- res$adjMaxP_pvals
reg_results_ILO$FDR_qval <- fdrq
View(reg_results_ILO)

reg_results_ILO$avg_Effect <- (reg_results_ILO$Est_Bracken + reg_results_ILO$Est_Metaphlan)/2

genus_s__species <- str_split_i(reg_results_ILO$tax.name_Bracken, "g__", 2)
reg_results_ILO$Species <- paste0(str_sub(genus_s__species, 1, 1), ". ", str_split_i(genus_s__species, "s__", 2))

# FDR Qvals for each method:
reg_results_ILO$FDR_Bracken <- p.adjust(reg_results_ILO$Pval_Bracken, method="BH")
reg_results_ILO$FDR_Metaphlan <- p.adjust(reg_results_ILO$Pval_Metaphlan, method="BH")

# add a column of NAs
reg_results_ILO$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_ILO$diffexpressed[reg_results_ILO$avg_Effect > 0 & reg_results_ILO$FDR_qval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_ILO$diffexpressed[reg_results_ILO$avg_Effect < 0 & reg_results_ILO$FDR_qval < 0.05] <- "DOWN"

reg_results_ILO$delabel <- NA
reg_results_ILO$delabel[reg_results_ILO$diffexpressed != "NO"] <- reg_results_ILO$Species[reg_results_ILO$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_ILO, aes(x=avg_Effect, y=-log10(FDR_qval), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
#  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Combined_volcano.png"))

# add a column of NAs
reg_results_ILO$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_ILO$diffexpressed[reg_results_ILO$Est_Bracken > 0 & reg_results_ILO$Pval_Bracken < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_ILO$diffexpressed[reg_results_ILO$Est_Bracken < 0 & reg_results_ILO$Pval_Bracken < 0.05] <- "DOWN"

reg_results_ILO$delabel <- NA
reg_results_ILO$delabel[reg_results_ILO$diffexpressed != "NO"] <- reg_results_ILO$Species[reg_results_ILO$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_ILO, aes(x=Est_Bracken, y=-log10(Pval_Bracken), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Bracken_volcano.png"))

# add a column of NAs
reg_results_ILO$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_ILO$diffexpressed[reg_results_ILO$Est_Metaphlan > 0 & reg_results_ILO$Pval_Metaphlan < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_ILO$diffexpressed[reg_results_ILO$Est_Metaphlan < 0 & reg_results_ILO$Pval_Metaphlan < 0.05] <- "DOWN"

reg_results_ILO$delabel <- NA
reg_results_ILO$delabel[reg_results_ILO$diffexpressed != "NO"] <- reg_results_ILO$Species[reg_results_ILO$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_ILO, aes(x=Est_Metaphlan, y=-log10(Pval_Metaphlan), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Metaphlan_volcano.png"))


## Any significant results with conflicting directions?
reg_results_ILO %>% filter(FDR_qval<0.05 & sign(Est_Bracken) != sign(Est_Metaphlan)) %>% pull(Species)


reg_results_ILO_long <- reg_results_ILO %>%
  pivot_longer(
    cols=c(Est_Bracken, Est_Metaphlan, Pval_Bracken, Pval_Metaphlan, FDR_Bracken, FDR_Metaphlan),
    names_to = c(".value", "Method"),
    names_pattern = "([A-Za-z]+)_([A-Za-z]+)"
  ) 


# plot adding up all layers we have seen so far
ggplot(data=reg_results_ILO_long, aes(x=Est, y=-log10(FDR), label=delabel, col=Method)) +
    geom_line(aes(group=tax.ID), col="gray") +
  geom_point(size=3, alpha=0.5) + 
  theme_classic() +
 # geom_text_repel() +
 # scale_color_manual(values=c("green", "purple")) +
  scale_color_viridis_d(begin = .2, end=.8)+
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")+ xlab("Estimate")
#ggsave(file.path(outdir, "BothMethods_volcano.png"), height=3, width=5, dpi=600)


ggplot(reg_results_ILO, aes(x=-log(FDR_Bracken), y=-log(FDR_Metaphlan), label=Species)) + geom_point() + geom_abline(slope=1) +
  geom_hline(yintercept = -log(0.05)) + geom_vline(xintercept = -log(0.05)) + theme_classic() + geom_text_repel() 
#ggsave(file.path(outdir, "FDR_bothMethods.png"))


ggplot(reg_results_ILO, aes(x=Est_Bracken, y=Est_Metaphlan, label=Species)) + geom_point() + geom_abline(slope=1)  + theme_classic() +
  geom_text_repel() 
#ggsave(file.path(outdir, "Estimates_bothMethods.png"))


### Xu et al Cohort ###


reg_results_Xu <- readxl::read_excel(file.path(indir, "Xu_Common_taxa_regression.xlsx"), sheet = "R1K3")
colnames(reg_results_Xu)

pvals <- as.matrix(reg_results_Xu %>% dplyr::select(Pval_Bracken, Pval_Metaphlan))
dim(pvals)

### Run AdjMaxP
res <- runAdjMaxP(pvals)

### AdjMaxP p-values
hist(res$adjMaxP_pvals, breaks=100)
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

# Investigate odd pattern
hist(reg_results_Xu$Pval_Bracken, breaks = 1000)
hist(reg_results_Xu$Pval_Metaphlan, breaks = 1000)

fdrq <- p.adjust(res$adjMaxP_pvals, method="BH")

reg_results_Xu$AdjMaxP <- res$adjMaxP_pvals
reg_results_Xu$FDR_qval <- fdrq
View(reg_results_Xu)

reg_results_Xu$avg_Effect <- (reg_results_Xu$Est_Bracken + reg_results_Xu$Est_Metaphlan)/2

genus_s__species <- str_split_i(reg_results_Xu$tax.name_Bracken, "g__", 2)
reg_results_Xu$Species <- paste0(str_sub(genus_s__species, 1, 1), ". ", str_split_i(genus_s__species, "s__", 2))

# FDR Qvals for each method:
reg_results_Xu$FDR_Bracken <- p.adjust(reg_results_Xu$Pval_Bracken, method="BH")
reg_results_Xu$FDR_Metaphlan <- p.adjust(reg_results_Xu$Pval_Metaphlan, method="BH")

# add a column of NAs
reg_results_Xu$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_Xu$diffexpressed[reg_results_Xu$avg_Effect > 0 & reg_results_Xu$FDR_qval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_Xu$diffexpressed[reg_results_Xu$avg_Effect < 0 & reg_results_Xu$FDR_qval < 0.05] <- "DOWN"

reg_results_Xu$delabel <- NA
reg_results_Xu$delabel[reg_results_Xu$diffexpressed != "NO"] <- reg_results_Xu$Species[reg_results_Xu$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_Xu, aes(x=avg_Effect, y=-log10(FDR_qval), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Combined_volcano_Xu.png"))

# add a column of NAs
reg_results_Xu$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_Xu$diffexpressed[reg_results_Xu$Est_Bracken > 0 & reg_results_Xu$Pval_Bracken < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_Xu$diffexpressed[reg_results_Xu$Est_Bracken < 0 & reg_results_Xu$Pval_Bracken < 0.05] <- "DOWN"

reg_results_Xu$delabel <- NA
reg_results_Xu$delabel[reg_results_Xu$diffexpressed != "NO"] <- reg_results_Xu$Species[reg_results_Xu$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_Xu, aes(x=Est_Bracken, y=-log10(Pval_Bracken), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Bracken_volcano_Xu.png"))

# add a column of NAs
reg_results_Xu$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_results_Xu$diffexpressed[reg_results_Xu$Est_Metaphlan > 0 & reg_results_Xu$Pval_Metaphlan < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_results_Xu$diffexpressed[reg_results_Xu$Est_Metaphlan < 0 & reg_results_Xu$Pval_Metaphlan < 0.05] <- "DOWN"

reg_results_Xu$delabel <- NA
reg_results_Xu$delabel[reg_results_Xu$diffexpressed != "NO"] <- reg_results_Xu$Species[reg_results_Xu$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
ggplot(data=reg_results_Xu, aes(x=Est_Metaphlan, y=-log10(Pval_Metaphlan), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey")
#ggsave(file.path(outdir, "Metaphlan_volcano_Xu.png"))


## Any significant results with conflicting directions? How does Eric's method handle this?
reg_results_Xu %>% filter(FDR_qval<0.05 & sign(Est_Bracken) != sign(Est_Metaphlan)) %>% pull(Species)
#character(0)


reg_results_Xu_long <- reg_results_Xu %>%
  pivot_longer(
    cols=c(Est_Bracken, Est_Metaphlan, Pval_Bracken, Pval_Metaphlan, FDR_Bracken, FDR_Metaphlan),
    names_to = c(".value", "Method"),
    names_pattern = "([A-Za-z]+)_([A-Za-z]+)"
  ) 




# plot adding up all layers we have seen so far
ggplot(data=reg_results_Xu_long, aes(x=Est, y=-log10(FDR), label=delabel, col=Method)) +
  geom_line(aes(group=tax.ID), col="gray") +
  geom_point(size=3, alpha=0.5) + 
  theme_classic() +
  # geom_text_repel() +
  # scale_color_manual(values=c("green", "purple")) +
  scale_color_viridis_d(begin = .2, end=.8)+
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey") +xlab("Estimate")
#ggsave(file.path(outdir, "BothMethods_volcano_Xu.png"), height=3, width=5, dpi=600)


ggplot(reg_results_Xu, aes(x=-log(FDR_Bracken), y=-log(FDR_Metaphlan), label=Species)) + geom_point() + geom_abline(slope=1) +
  geom_hline(yintercept = -log(0.05)) + geom_vline(xintercept = -log(0.05)) + theme_classic() + geom_text_repel()
#ggsave(file.path(outdir, "FDR_bothMethods_Xu.png"))


ggplot(reg_results_Xu, aes(x=Est_Bracken, y=Est_Metaphlan, label=Species)) + geom_point() + geom_abline(slope=1)  + theme_classic() +
  geom_text_repel() 
#ggsave(file.path(outdir, "Estimates_bothMethods_Xu.png"))


### Read in full ILO results:

reg_all_ILO_bkn <- readRDS(file.path(indir, "ILO_BrackenBracken_parsed_results.rds"))[["species"]][[9]][["data"]]
reg_all_ILO_mp4 <- readRDS(file.path(indir, "ILO_MetaphlanMetaphlan_parsed_results.rds"))[["species"]]

reg_results_pvals <- reg_results_ILO %>% dplyr::select(tax.ID, AdjMaxP)

reg_union <- full_join(reg_all_ILO_bkn, reg_all_ILO_mp4, by="tax.ID", suffix=c(".Bracken", ".Metaphlan"))
reg_union <- left_join(reg_union, reg_results_pvals)
reg_union <- reg_union %>% 
  mutate(
    Method = case_when(
      is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ "Metaphlan only",
      !is.na(Est.Bracken) & is.na(Est.Metaphlan) ~ "Bracken only",
      !is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ "Both",
      .default = "Other"
    ),
    Pval_combined = case_when(
      is.na(Pval.Bracken) & !is.na(Pval.Metaphlan) ~ Pval.Metaphlan,
      !is.na(Pval.Bracken) & is.na(Pval.Metaphlan) ~ Pval.Bracken,
      !is.na(Pval.Bracken) & !is.na(Pval.Metaphlan) ~ AdjMaxP,
      .default = NA
    ),
    Taxname_combined = case_when(
      !is.na(Est.Bracken) ~ tax.name.Bracken,
      is.na(Est.Bracken)  ~ tax.name.Metaphlan,
      .default = ""
    ),
    Est_combined = case_when(
      is.na(Est.Bracken)  ~ Est.Metaphlan,
      !is.na(Est.Bracken) & is.na(Est.Metaphlan) ~ Est.Bracken,
      !is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ (Est.Metaphlan + Est.Bracken)/2,
    )
  )
reg_union$FDR_combined <- p.adjust(reg_union$Pval_combined, method="BH")

reg_union$genus_s__species <- NA
reg_union$Species <- NA
reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)] <- str_split_i(reg_union$Taxname_combined[!is.na(reg_union$Est.Bracken)], "g__", 2)
reg_union$Species[!is.na(reg_union$Est.Bracken)] <- paste0(str_split_i(reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)], "s__", 1), str_split_i(reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)], "s__", 2))
reg_union$Species[is.na(reg_union$Est.Bracken)] <- reg_union$taxaname.Metaphlan[is.na(reg_union$Est.Bracken)]

# Save for different upset plot:
reg_union2 <- reg_union


# add a column of NAs
reg_union$diffexpressed <- "Not significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05] <- "Up in Age (Based on both methods)"
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Bracken)] <- "Up in Age (Based on Metaphlan only)"
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Metaphlan)] <- "Up in Age (Based on Bracken only)"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05] <- "Down in Age (Based on both methods)"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Bracken)] <- "Down in Age (Based on Metaphlan only)"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Metaphlan)] <- "Down in Age (Based on Bracken only)"

reg_union_ILO <- reg_union



### Read in full Xu results:

reg_all_Xu_bkn <- readRDS(file.path(indir, "Asian_BrackenBracken_parsed_results.rds"))[["species"]][[6]][["data"]]
reg_all_Xu_mp4 <- readRDS(file.path(indir, "Asian_MetaphlanMetaphlan_parsed_results.rds"))[["species"]]

reg_results_pvals <- reg_results_Xu %>% dplyr::select(tax.ID, AdjMaxP)

reg_union <- full_join(reg_all_Xu_bkn, reg_all_Xu_mp4, by="tax.ID", suffix=c(".Bracken", ".Metaphlan"))
reg_union <- left_join(reg_union, reg_results_pvals)
reg_union <- reg_union %>% 
  mutate(
    Method = case_when(
      is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ "Metaphlan only",
      !is.na(Est.Bracken) & is.na(Est.Metaphlan) ~ "Bracken only",
      !is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ "Both",
      .default = "Other"
    ),
    Pval_combined = case_when(
      is.na(Pval.Bracken) & !is.na(Pval.Metaphlan) ~ Pval.Metaphlan,
      !is.na(Pval.Bracken) & is.na(Pval.Metaphlan) ~ Pval.Bracken,
      !is.na(Pval.Bracken) & !is.na(Pval.Metaphlan) ~ AdjMaxP,
      .default = NA
    ),
    Taxname_combined = case_when(
      !is.na(Est.Bracken) ~ tax.name.Bracken,
      is.na(Est.Bracken)  ~ tax.name.Metaphlan,
      .default = ""
    ),
    Est_combined = case_when(
      is.na(Est.Bracken)  ~ Est.Metaphlan,
      !is.na(Est.Bracken) & is.na(Est.Metaphlan) ~ Est.Bracken,
      !is.na(Est.Bracken) & !is.na(Est.Metaphlan) ~ (Est.Metaphlan + Est.Bracken)/2,
    )
  )
reg_union$FDR_combined <- p.adjust(reg_union$Pval_combined, method="BH")

reg_union$genus_s__species <- NA
reg_union$Species <- NA
reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)] <- str_split_i(reg_union$Taxname_combined[!is.na(reg_union$Est.Bracken)], "g__", 2)
reg_union$Species[!is.na(reg_union$Est.Bracken)] <- paste0(str_split_i(reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)], "s__", 1), str_split_i(reg_union$genus_s__species[!is.na(reg_union$Est.Bracken)], "s__", 2))
reg_union$Species[is.na(reg_union$Est.Bracken)] <- reg_union$taxaname.Metaphlan[is.na(reg_union$Est.Bracken)]


# add a column of NAs
reg_union$diffexpressed <- "Not significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05] <- "Up in Age (Based on both methods)"
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Bracken)] <- "Up in Age (Based on Metaphlan only)"
reg_union$diffexpressed[reg_union$Est_combined > 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Metaphlan)] <- "Up in Age (Based on Bracken only)"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05] <- "Down in Age (Based on both methods)"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Bracken)] <- "Down in Age (Based on Metaphlan only)"
reg_union$diffexpressed[reg_union$Est_combined < 0 & reg_union$FDR_combined < 0.05 & is.na(reg_union$Est.Metaphlan)] <- "Down in Age (Based on Bracken only)"

reg_union_Xu <- reg_union




### Plot both cohorts' volcano plots
reg_union <- reg_union_ILO

reg_union$delabel <- NA
reg_union$delabel[which(reg_union$FDR_combined <0.05 & reg_union$tax.ID %in% reg_union_Xu$tax.ID[reg_union_Xu$FDR_combined<0.05])] <-
  reg_union$Species[which(reg_union$FDR_combined <0.05 & reg_union$tax.ID %in% reg_union_Xu$tax.ID[reg_union_Xu$FDR_combined<0.05])]
reg_union$delabel[reg_union$FDR_combined<10^(-4)] <- reg_union$Species[reg_union$FDR_combined<10^(-4)]

ggplot(data=reg_union, aes(x=Est_combined, y=-log10(FDR_combined), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() +
  xlim(-0.12, 0.12) +
  geom_text_repel(max.overlaps = 7, min.segment.length = 0.001, size = 3.5) +
  scale_color_manual(values=c("blue", "lightblue", "lightgreen", "darkgrey", "red", "pink2", "orange")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.05), col="darkgrey") + labs(col = "Differential Expression by Age \nvia Correlated Meta-Analysis") + xlab("Combined Estimate") + ylab("-log10(combined FDR)")
ggsave(file.path(outdir, "ILO_volcano_union_20250214.png"), height=4, width=9, dpi=600)



reg_union <- reg_union_Xu

reg_union$delabel <- NA
reg_union$delabel[which(reg_union$FDR_combined <0.05 & reg_union$tax.ID %in% reg_union_ILO$tax.ID[reg_union_ILO$FDR_combined<0.05])] <-
  reg_union$Species[which(reg_union$FDR_combined <0.05 & reg_union$tax.ID %in% reg_union_ILO$tax.ID[reg_union_ILO$FDR_combined<0.05])]
reg_union$delabel[reg_union$FDR_combined<10^(-5)] <- reg_union$Species[reg_union$FDR_combined<10^(-5)]


vh_line <- data.frame(
  yintercept = -log10(0.05),
  FDR_combined = 12
)


# plot adding up all layers we have seen so far
ggplot(data=reg_union, aes(x=Est_combined, y=-log10(FDR_combined), label=delabel, col=diffexpressed)) +
  geom_point() + 
  theme_classic() + 
  xlim(-0.12, 0.12) +
  geom_text_repel(max.overlaps = 20, min.segment.length = 0.001, size = 3.5) +
  scale_color_manual(values=c("blue", "lightblue", "lightgreen", "darkgrey", "red", "pink2", "orange")) +
  #  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  labs(col = "Differential Expression by Age \nvia Correlated Meta-Analysis") + 
  xlab("Combined Estimate") + ylab("-log10(combined FDR)") +

#geom_blank(data = data.frame(x = Inf, y = c(6, 10))) +
#  geom_text_repel(aes(label = label)) +
  facet_grid(-log10(FDR_combined) <= 6 ~ ., space = "free_y", scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
geom_hline(data = vh_line, aes(yintercept = yintercept), col="darkgrey") 


ggsave(file.path(outdir, "Xu_volcano_union_20250214.png"), height=4, width=9, dpi=600)

##### Combine data for upset plots

union_bothcohorts <- full_join(reg_union_ILO, reg_union_Xu, by=c("tax.ID", "Species"), suffix=c(".ILO", ".Xu"))
# Save for local plotting
write.csv(union_bothcohorts, "volcano_plots/union_bothcohorts.csv")


View(union_bothcohorts %>% dplyr::select(tax.ID, Species, diffexpressed.ILO, diffexpressed.Xu))
View(union_bothcohorts %>% dplyr::filter(diffexpressed.ILO==diffexpressed.Xu & diffexpressed.ILO!="NO"))

table(union_bothcohorts$diffexpressed.ILO, union_bothcohorts$diffexpressed.Xu, useNA = "ifany")

## Simplify for tables to save
reg_union_ILO_tosave <- reg_union_ILO %>% 
  dplyr::select(c(tax.ID, Species,  tax.name.Bracken , Est.Bracken,  Pval.Bracken,
                                                 tax.name.Metaphlan, Est.Metaphlan, Pval.Metaphlan,   AdjMaxP, Method,
                                                 Est_combined, Pval_combined, FDR_combined)) %>%
  arrange(FDR_combined)
write.csv(reg_union_ILO_tosave, "volcano_plots/ILO_species_results.csv", row.names = F)

reg_union_Xu_tosave <- reg_union_Xu %>% 
  dplyr::select(c(tax.ID, Species,  tax.name.Bracken , Est.Bracken,  Pval.Bracken, 
                  tax.name.Metaphlan, Est.Metaphlan, Pval.Metaphlan,   AdjMaxP, Method,
                  Est_combined, Pval_combined, FDR_combined)) %>%
  arrange(FDR_combined)
write.csv(reg_union_Xu_tosave, "volcano_plots/Xu_species_results.csv", row.names = F)


