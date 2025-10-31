#load libraries
library(BiocManager)
library(dplyr)   
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(microbiome)
library(DESeq2) # differential abundance analysis
library(rstatix) # dunn_test
library(ANCOMBC)
library(DT)
library(ComplexHeatmap)# For advanced heatmaps
library(foreach)
library(rngtools)

print(paste("Libraries are loaded"))
args <- commandArgs(trailingOnly=TRUE)
comp <- c("mucus","tissue","skeleton")

print(paste("Looping over mucus tissue and skeleton compartments..."))
for (i in comp){
sink_name = paste0(i,".","functional_ancombc_analysis_log.txt")
sink(sink_name,append=FALSE,split=TRUE)

#Import from .qza file into a phyloseq object
funk_path_name <- paste0("./",i,"_picrust2_out_pipeline","/",i,"_picrust_pathway_table.qza")
asv <- qza_to_phyloseq(features = funk_path_name)

#### Import Metadata read.table
map_path_name <- paste0("./",i,"/",i,"_combined_metadata.txt")
metadata <- read.table(file =map_path_name ,header=T,comment.char="", sep="\t")

### Import taxonomy from biom output as .tsv format using read.table
taxonomy <- read.table(file = "metacyc_prokaryote_descrip.tsv", sep = "\t", header = T,row.names = 1)
taxonomy$Family <-row.names(taxonomy)
colnames(taxonomy)[1] <-"Phylum"
# create matrix format for OTU and taxonomy table
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(taxonomy))

# Set metadata
SAMPLE <- sample_data(metadata)
sample_names(SAMPLE) <- sample_names(OTU)

# Create Working phyloseq object
funk_phylo <- phyloseq(OTU,tax1,SAMPLE)
funk_phylo
## Assign sample data to meta_data variable
meta_data <- sample_data(funk_phylo)
# factor the neutrality to set neutral first
meta_data$neutrality = factor(meta_data$neutrality, levels = c("neutral","below","above"))
# reassign factored metadata to phyloseq object
sample_data(funk_phylo) = meta_data

## Subset between two groups
#subset_phylo <- subset_samples(funk_phylo, neutrality!=("neutral"))
#subset_phylo

######################## Ancombc analysis ##############################################
# comparison between Neutral and non-neutral taxa
print(paste("Running ancombc2 analysis on:",i))
set.seed(111)
output = ancombc2(data = funk_phylo, tax_level = "Phylum",
                  fix_formula = "neutrality", p_adj_method = "holm", pseudo_sens = TRUE,
                  group = "neutrality", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05,pairwise = TRUE)
## Pull out results
res = output$res_pair


df_fig_neutral = res %>%
  dplyr::filter(diff_neutralityabove == 1 | 
                  diff_neutralitybelow == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_neutralityabove == 1, 
                              round(lfc_neutralityabove, 2), 0),
                lfc2 = ifelse(diff_neutralitybelow == 1, 
                              round(lfc_neutralitybelow, 2), 0),
                lfc3 = ifelse(diff_neutralityabove_neutralitybelow == 1, 
                              round(lfc_neutralityabove_neutralitybelow, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) 


## Assign variable names 
df_fig_neutral$group = recode(df_fig_neutral$group, 
                              `lfc1` = "Above - Neutral",
                              `lfc2` = "Below - Neutral",
                              `lfc3` = "Above- Below")
df_fig_neutral$group = factor(df_fig_neutral$group, 
                              levels = c("Above - Neutral",
                                         "Below - Neutral",
                                         "Above- Below"))

## set figure range 
lo = floor(min(df_fig_neutral$value))
up = ceiling(max(df_fig_neutral$value))
mid = (lo + up)/2

fig_title = print(paste(i,"","pairwise log-fold change in predicted pathways between \n non-neutral and neutral communities"))

## Build Heatmap
fig_neutral = df_fig_neutral %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#006CD1", high = "#994F00", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = "black"), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = fig_title) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(paste("Creating Heatmap pdf File"))

## Print Neutral Model Plot
heat_name <- paste0(i,"_ancombc_heatmap.pdf")
ggsave(fig_neutral, filename=heat_name, width = 20, height = 10, dpi= 600)

}

print("Finished!")
