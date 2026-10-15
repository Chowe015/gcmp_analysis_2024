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
library(gplots)
library(vegan)
library(RVAideMemoire)
library(multcompView)

## Step 2 Run Sloan Neutral Model to produce neutral table and figures
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
setwd("C:/Users/Colin Howe/Desktop/gcmp_final_analysis")

## Comparison between above neutral and below categories for mucus tissue and skeleton
## Import mapping and otu functional table for each compartment and create a specific phyloseq object
## Mucus functional phyloseq

sink_name = paste0("Picrust2_ancombc_analysis_log.txt")
sink(sink_name,append=FALSE,split=TRUE)

mucus_asv <- read.table(file = "./final_picrust2/mucus_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_mucus <- mucus_asv[, colSums(mucus_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
mucus_map <- read.table(file ="./final_picrust2/mucus_picrust2_out_pipeline/mucus_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
muc_otu <- otu_table(as.matrix(clean_mucus), taxa_are_rows = TRUE)

# Set metadata
muc_sample <- sample_data(mucus_map)

# Create Working phyloseq object
muc_funk_ph <- phyloseq(muc_otu,muc_sample)
muc_funk_ph

## Tissue functional phyloseq
tissue_asv <- read.table(file = "./final_picrust2/tissue_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_tissue <- tissue_asv[, colSums(tissue_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
tissue_map <- read.table(file ="./final_picrust2/tissue_picrust2_out_pipeline/tissue_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
tis_otu <- otu_table(as.matrix(clean_tissue), taxa_are_rows = TRUE)

# Set metadata
tis_sample <- sample_data(tissue_map)


# Create Working phyloseq object
tis_funk_ph <- phyloseq(tis_otu,tis_sample)
tis_funk_ph

## skeleton functional phyloseq
skeleton_asv <- read.table(file = "./final_picrust2/skeleton_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_skeleton <- skeleton_asv[, colSums(skeleton_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
skeleton_map <- read.table(file ="./final_picrust2/skeleton_picrust2_out_pipeline/skeleton_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
ske_otu <- otu_table(as.matrix(clean_skeleton), taxa_are_rows = TRUE)

# Set metadata
ske_sample <- sample_data(skeleton_map)


# Create Working phyloseq object
ske_funk_ph <- phyloseq(ske_otu,ske_sample)
ske_funk_ph

## join each phyloseq object
merge_phyloseq(muc_funk_ph,tis_funk_ph,ske_funk_ph) -> merged_funk_ph

# Calculate DISTANCE MATRICES
## All samples
bray_merge = phyloseq::distance(merged_funk_ph, method="bray")

# Create data frame for downstream analysis
all_merge_df <- data.frame(sample_data(merged_funk_ph))

# join metacyc pathway descriptions as taxonomy 
taxonomy <- read.table(file = "./final_picrust2/metacyc_pathways_decrip_full.txt", sep = "\t", header = T,row.names = 1)
names(taxonomy)[1] <-paste("Phylum")
names(taxonomy)[2] <-paste("Class")
names(taxonomy)[3] <-paste("Order")

# Transform into a matrix
tax1 = tax_table(as.matrix(taxonomy))
tax_table(merged_funk_ph) <- tax1

merged_funk_ph

######################## Ancombc analysis ##############################################
# comparison between Neutral and non-neutral taxa
print(paste("Running ancombc2 analysis on Merged Phyloseq Object"))
set.seed(111)

# Set groups are factor
sample_data(merged_funk_ph)$compartment_neutrality <- as.factor(sample_data(merged_funk_ph)$compartment_neutrality)

output = ancombc2(data = merged_funk_ph, tax_level = "Class", fix_formula = "compartment_neutrality", p_adj_method = "fdr", pseudo_sens = TRUE, group = "compartment_neutrality", global=TRUE, pairwise=TRUE, struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05)

## Pull out results
res_global = output$res_global

res = output$res
res_pairwise = output$res_pair
colnames(res)[2] <- "lfc_Intercept"
colnames(res)[38] <- "q_Intercept"
write.table(res, file="./final_picrust2/ancombc2_picrust2_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)
write.table(res_pairwise, file="./final_picrust2/ancombc2_picrust2_pairwise_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)
write.table(res_global, file="./final_picrust2/ancombc2_picrust2_global_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)


print(paste("Organizing results and generating Bubble Plot"))
df_fig_neutral = res %>%
  dplyr::filter(lfc_Intercept == 1 |diff_compartment_neutralitymucus_below == 1 |diff_compartment_neutralitymucus_neutral == 1 |diff_compartment_neutralitytissue_above == 1 |diff_compartment_neutralitytissue_below == 1 |diff_compartment_neutralitytissue_neutral == 1 |diff_compartment_neutralityskeleton_above == 1 |diff_compartment_neutralityskeleton_below == 1 |diff_compartment_neutralityskeleton_neutral == 1) %>%
  dplyr::mutate(lfc1 = ifelse(q_Intercept == 1, 
                              round(lfc_Intercept, 2), 0),
                lfc2 = ifelse(diff_compartment_neutralitymucus_below == 1, 
                              round(lfc_compartment_neutralitymucus_below, 2), 0),
                lfc3 = ifelse(diff_compartment_neutralitymucus_neutral == 1, 
                              round(lfc_compartment_neutralitymucus_neutral, 2), 0),
                lfc4 = ifelse(diff_compartment_neutralitytissue_above == 1, 
                              round(lfc_compartment_neutralitytissue_above, 2), 0),
                lfc5 = ifelse(diff_compartment_neutralitymucus_below == 1, 
                              round(lfc_compartment_neutralitytissue_below, 2), 0),
                lfc6 = ifelse(diff_compartment_neutralitytissue_neutral == 1, 
                              round(lfc_compartment_neutralitytissue_neutral, 2), 0),
                lfc7 = ifelse(diff_compartment_neutralityskeleton_above == 1, 
                              round(lfc_compartment_neutralityskeleton_above, 2), 0),
                lfc8 = ifelse(diff_compartment_neutralityskeleton_below == 1, 
                              round(lfc_compartment_neutralityskeleton_below, 2), 0),
                lfc9 = ifelse(diff_compartment_neutralityskeleton_neutral == 1, 
                              round(lfc_compartment_neutralityskeleton_neutral, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc9, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_neutral$group = recode(df_fig_neutral$group, 
                              `lfc1` = "mucus_above",
                              `lfc2` = "mucus_below",
                              `lfc3` = "mucus_neutral",
                              `lfc4` = "tissue_above",
                              `lfc5` = "tissue_below",
                              `lfc6` = "tissue_neutral",
                              `lfc7` = "skeleton_above",
                              `lfc8` = "skeleton_below",
                              `lfc9` = "skeleton_neutral")

df_fig_neutral$group = factor(df_fig_neutral$group, 
                              levels=c("mucus_above",
                                       "mucus_neutral",
                                       "mucus_below",
                                       "tissue_above",
                                       "tissue_neutral",
                                       "tissue_below",
                                       "skeleton_above",
                                       "skeleton_neutral",
                                       "skeleton_below"))


#### Bubble Figure
fig_neutral = df_fig_neutral %>%
  ggplot(aes(x = group, y = taxon)) + theme_minimal() +
  geom_vline(xintercept=0, linetype= "dashed", color = "grey50", linewidth= 0.5) +
  geom_point(aes(size=abs(value), fill= value >=0), shape=21, color="black", alpha=0.8) +
    scale_size_continuous(
    range = c(1, 10), 
    name = "Absolute Value"
  ) + scale_fill_manual(
    values = c("TRUE" = "#006CD1", "FALSE" = "#994F00"),
    labels = c("FALSE" = "Negative Lfc", "TRUE" = "Positive Lfc"),
    name = "Logfold Change")

print(paste("Printing Ancombc2 Analysis Results &Bubble Plot"))
#write.table(df_fig_neutral, file="./final_picrust2/ancombc2_bubble_plot_results.tsv", row.names=TRUE, sep = "\t")
ggsave("./final_picrust2/Picrust2 Compartment-Neutral Bubble Plot.jpg", fig_neutral, width = 20, height = 10, dpi= 600)


##### Multicomp tool 


# Identify columns containing the adjusted p-values (q-values)
# ANCOM-BC2 names these as 'q_groupLevel1_groupLevel2'
q_cols <- grep("^q_", colnames(res_pair), value = TRUE)

# --- Step 2: Function to convert pairwise results to letters for a single taxon ---
get_taxon_letters <- function(taxon_row, q_columns) {
  # Extract q-values for this specific taxon row
  p_vals <- as.numeric(taxon_row[q_columns])
  
  # Clean up column names to get exact pairwise group strings (e.g., "groupA-groupB")
  # This strips the leading "q_" or "q_mainvar" added by ANCOMBC
  names(p_vals) <- gsub("^q_", "", q_columns) 
  
  # Convert the named p-value vector into a significance group letter object
  # If any p-value is NA (due to structural zeros), we handle it using tryCatch
  letters_obj <- tryCatch({
    multcompLetters(p_vals, threshold = 0.05)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(letters_obj)) return(NULL)
  
  # Convert to a clean data frame: Groups and their assigned Letters
  df_letters <- data.frame(
    Group = names(letters_obj$Letters),
    Letter = letters_obj$Letters,
    stringsAsFactors = FALSE
  )
  return(df_letters)
}

# --- Step 3: Loop through all taxa and build the letters database ---
letters_list <- list()

for (i in 1:nrow(res_pair)) {
  taxon_name <- res_pair$taxon[i]
  taxon_letters <- get_taxon_letters(res_pair[i, ], q_cols)
  
  if (!is.null(taxon_letters)) {
    taxon_letters$taxon <- taxon_name
    letters_list[[i]] <- taxon_letters
  }
}

# Combine all rows into a single long-format data frame
df_all_letters <- bind_rows(letters_list)

#################### Analysis of Variance (ANOVA) #####################################
print(paste0("Compartment and Neutrality anova analysis"))
set.seed(129)

adonis2(bray_merge~compartment_neutrality, data=all_merge_df) -> result_all
print(result_all)

print(paste0("Printing Anova and permanova Result"))
bray_comp_pairwise <-pairwise.perm.manova(bray_merge, all_merge_df$compartment_neutrality,nperm=1000) 
print(bray_comp_pairwise)

#Analyze multivariate homogeneity of group dispersions (centroids)
dispersion <- betadisper(bray_merge, group = all_merge_df$compartment_neutrality)

# 7. Perform Tukey's HSD on the dispersion object
tukey_output <- TukeyHSD(dispersion)
as.data.frame(tukey_output$group)-> tukey_result
# 8. View the pairwise comparison results
print(tukey_output)

# 9. Plot the Tukey HSD results
plot(tukey_output)

write.table(tukey_result, file="./final_picrust2/merged_comp_neutrality_results.tsv", row.names = TRUE,sep="\t")


####################### Anova Comparison ####################################
print(paste("Running Anova for Picrust2 functional results for Compartment and Neutrality"))
comp <- c("mucus","tissue","skeleton")
#i = "skeleton"

for (i in comp){ 
  
  funk_path_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv")
  #asv <- qza_to_phyloseq(features = funk_path_name)
  asv_test <- read.table(file = funk_path_name,header=T,comment.char="", sep="\t", row.names = 1)
  df_clean <- asv_test[, colSums(asv_test != 0) != 0, drop = FALSE]
  
  #### Import Metadata read.table
  map_path_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline","/",i,"_combined_metadata.txt")
  metadata <- read.table(file =map_path_name ,header=T,comment.char="", sep="\t",row.names = 1)
  
  ### Import taxonomy from biom output as .tsv format using read.table
  taxonomy <- read.table(file = "./final_picrust2/metacyc_pathways_decrip_full.txt", sep = "\t", header = T,row.names = 1)
  names(taxonomy)[1] <-paste("Phylum")
  names(taxonomy)[2] <-paste("Class")
  names(taxonomy)[3] <-paste("Order")
  
  # create matrix format for OTU and taxonomy table
  OTU <- otu_table(as.matrix(df_clean), taxa_are_rows = TRUE)
  tax1 = tax_table(as.matrix(taxonomy))
  
  # Set metadata
  SAMPLE <- sample_data(metadata)
  
  
  # Create Working phyloseq object
  funk_phylo <- phyloseq(OTU,tax1,SAMPLE)
  funk_phylo
  
  
  # Calculate DISTANCE MATRICES
  ## All samples
  bray_all = phyloseq::distance(funk_phylo, method="bray")
  
  # Create data frame for downstream analysis
  all_df <- data.frame(sample_data(funk_phylo))
  
  print(paste0(i," anova analysis"))
  set.seed(129)
  adonis2(bray_all~neutrality, data=all_df) -> result
  print(result)
  
  print(paste0(i," anova and permanova analysis"))
  bray_comp_pairwise <-pairwise.perm.manova(bray_all, all_df$neutrality,nperm=1000) 
  print(bray_comp_pairwise)
  
  #Analyze multivariate homogeneity of group dispersions (centroids)
  dispersion <- betadisper(bray_all, group = all_df$neutrality)
  
  # 7. Perform Tukey's HSD on the dispersion object
  tukey_output <- TukeyHSD(dispersion)
  
  # 8. View the pairwise comparison results
  print(tukey_output)
  
  # 9. Plot the Tukey HSD results
  plot(tukey_output)
  
}