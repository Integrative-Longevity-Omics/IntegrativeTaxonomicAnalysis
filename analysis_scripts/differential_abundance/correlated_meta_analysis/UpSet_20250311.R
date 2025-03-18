library(ComplexUpset)
library(ggplot2)
library(dplyr)

outdir <- "/restricted/projectnb/uh2-sebas/analysis/metagenomics/meg_analyses/volcano_plots"


union_bothcohorts <- read.csv(file.path(outdir, "union_bothcohorts.csv"), row.names = 1)


union_small <- union_bothcohorts %>% 
  dplyr::filter(p.adj.Metaphlan.ILO<0.05 | p.adj.Bracken.ILO <0.05 | p.adj.Metaphlan.Xu<0.05 | p.adj.Bracken.Xu <0.05) %>% 
  select(tax.ID, p.adj.Metaphlan.ILO, p.adj.Bracken.ILO, p.adj.Metaphlan.Xu, p.adj.Bracken.Xu,
         FDR_combined.ILO, FDR_combined.Xu)

union_small$ILO.Metaphlan <- union_small$p.adj.Metaphlan.ILO<0.05
union_small$ILO.Bracken <- union_small$p.adj.Bracken.ILO<0.05
union_small$Xu.Metaphlan <- union_small$p.adj.Metaphlan.Xu<0.05 
union_small$Xu.Bracken <- union_small$p.adj.Bracken.Xu<0.05 

union_small <- union_small %>%
  mutate(combined_cat = case_when(
    is.na(FDR_combined.ILO) | is.na(FDR_combined.Xu) ~ "Not identified by both methods",
    FDR_combined.ILO<0.05 & FDR_combined.Xu < 0.05 ~ "Both Cohorts",
    FDR_combined.ILO<0.05  ~ "ILO",
    FDR_combined.Xu < 0.05 ~ "Xu",
    .default = "Neither"
  ))

union_small$combined_cat.f <- factor(union_small$combined_cat, levels=c("Not identified by both methods", "Neither", "Both Cohorts", "ILO", "Xu" ))

vars1 <- c("ILO.Metaphlan", "ILO.Bracken", "Xu.Metaphlan", "Xu.Bracken")


union_small_Xu <- union_bothcohorts %>% filter(FDR_combined.Xu<0.05 | p.adj.Metaphlan.Xu<0.05 | p.adj.Bracken.Xu<0.05)
union_small_Xu$Significant.Metaphlan <- union_small_Xu$p.adj.Metaphlan.Xu<0.05 
union_small_Xu$Significant.Bracken <- union_small_Xu$p.adj.Bracken.Xu<0.05 

union_small_Xu <- union_small_Xu %>% 
  mutate(combined_cat = case_when(
    #is.na(p.adj.Metaphlan.Xu) | is.na(p.adj.Bracken.Xu) ~ "Not identified by both methods",
    FDR_combined.Xu<0.05 ~ "Significant (n=71)",
    FDR_combined.Xu>=0.05 ~ "Not Significant (n=18)"
  ))

p3 <- ComplexUpset::upset( 
  as.data.frame(union_small_Xu),
  c("Significant.Metaphlan", "Significant.Bracken"),
  name="",
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))# + scale_y_continuous(breaks=c(5,7.5, 20, 25))
  ),
  
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=T,
      mapping=aes(fill=combined_cat)
    )  + labs(fill = "Correlated \nMeta-analysis") 
  )
)
p3
png(file=file.path(outdir, "UpsetPlot_Xu_20250213.png"), width = 6, height = 3, res = 600, units = "in") 
p3
dev.off()

write.csv(union_small_Xu, file = file.path(outdir, "Xu_data_forUpsetPlot.csv"), row.names=F)

union_small_ILO <- union_bothcohorts %>% filter(FDR_combined.ILO<0.05 | p.adj.Metaphlan.ILO<0.05 | p.adj.Bracken.ILO<0.05)
union_small_ILO$Significant.Metaphlan <- union_small_ILO$p.adj.Metaphlan.ILO<0.05 
union_small_ILO$Significant.Bracken <- union_small_ILO$p.adj.Bracken.ILO<0.05 

union_small_ILO <- union_small_ILO %>% 
  mutate(combined_cat = case_when(
    #is.na(p.adj.Metaphlan.ILO) | is.na(p.adj.Bracken.ILO) ~ "Not identified by both methods",
    FDR_combined.ILO<0.05 ~ "Significant (n=87)",
    FDR_combined.ILO>=0.05 ~ "Not Significant (n=24)"
  ))

p4 <- ComplexUpset::upset( 
  as.data.frame(union_small_ILO),
  c("Significant.Metaphlan", "Significant.Bracken"),
  name="",
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))# + scale_y_continuous(breaks=c(5,7.5, 20, 25))
  ),
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=T,
      mapping=aes(fill=combined_cat)
    )  + labs(fill = "Correlated \nMeta-analysis") 
  )
)
p4
png(file=file.path(outdir, "UpsetPlot_ILO_20250213.png"), width = 6, height = 3, res = 600, units = "in") 
p4
dev.off()


write.csv(union_small_ILO, file = file.path(outdir, "ILO_data_forUpsetPlot.csv"), row.names = F)


library(ggvenn)

x1 <- list("ILO" = union_bothcohorts$tax.ID[which(union_bothcohorts$FDR_combined.ILO<0.05)],
     "Xu et al." = union_bothcohorts$tax.ID[which(union_bothcohorts$FDR_combined.Xu<0.05)])

p5 <- ggvenn(
  x1, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4.5, text_size = 4.5, 
)
p5
png(file=file.path(outdir, "Venn_20250213.png"), width = 4, height = 3, res = 600, units = "in") 
p5
dev.off()

# Table of species
union_combined <- union_bothcohorts %>% 
  filter(FDR_combined.ILO<0.05 | FDR_combined.Xu<0.05) %>%
  select(Species, tax.ID, 
         Est.Metaphlan.ILO, Pval.Metaphlan.ILO, p.adj.Metaphlan.ILO, 
         Est.Bracken.ILO, Pval.Bracken.ILO, p.adj.Bracken.ILO, 
         AdjMaxP.ILO, Pval_combined.ILO, FDR_combined.ILO,
         Est.Metaphlan.Xu, Pval.Metaphlan.Xu, p.adj.Metaphlan.Xu, 
         Est.Bracken.Xu, Pval.Bracken.Xu, p.adj.Bracken.Xu, 
         AdjMaxP.Xu, Pval_combined.Xu, FDR_combined.Xu)



union_17 <- union_combined %>% filter(FDR_combined.ILO<0.05 & FDR_combined.Xu<0.05) %>% arrange(FDR_combined.ILO) 

write.csv(union_17, file = file.path(outdir, "replicated_sp_acrossCohorts.csv"))

