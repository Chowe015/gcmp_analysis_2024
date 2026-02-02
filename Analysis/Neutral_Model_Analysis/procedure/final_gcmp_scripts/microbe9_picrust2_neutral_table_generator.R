#Loading Required libraries
library(tidyverse)
library(phyloseq)
library(BiocManager)
library(dplyr)
library(data.table)
library(readr)
library(ggpicrust2)
library(tibble)
library(ggprism)
library(patchwork)


## STep 6 Import picrust2 results from full analysis of gcmp samples and merge table to understand functional annotation
sink("Picrust2_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

print(paste("importing data..."))
neutral_results_path <-args[1]
kegg_path <-args[2]
pathway_path <-args[3]
biosample <-args[4]

#### Read neutral table input
neutral_results <- read.table(neutral_results_path, header=T, sep="\t", check.names = FALSE)
#neutral_results <- read.table("./Mucus/Mucus_neutral_model.tsv", header=T, sep="\t", check.names = FALSE)

## Write neutral tax table
#neutral_tax_table_name <- paste0(biosample,"all_neutral_tax_table.tsv")
#write.table(neutral_tax_table, file =neutral_tax_table,sep="\t",row.names = FALSE, col.names=TRUE)

## Filter for significant families
print(paste("Filtering Significant ASV from Neutral Table"))
#filter(neutral_results, padj <= 0.05) %>% select(c(id,padj,Phylum,Class,Order,Family,model)) -> sig_neutral_table

print(paste("Importing KO and Kegg contribution tables"))
pathway <- read.table(gzfile(pathway_path), header =T, sep = "\t")#,check.names = FALSE)
colnames(pathway)[which(names(pathway)== "taxon")] <- "id"

#ko_function <- read.table(gzfile(kegg_path), header =T, sep = "\t")#, check.names = FALSE)
#colnames(ko_function)[which(names(ko_function)== "taxon")] <- "id"

#### Inner_join feature table with significant microbial families 
print(paste("joining neutral table and ref taxonomy"))
non_neutral_path_contrib <- inner_join(pathway,neutral_results, by = "id")

neutral_pathway_name <- paste0(biosample,"_neutral_pathway_contrib.tsv.gz")
write.table(non_neutral_path_contrib, gzfile(neutral_pathway_name),sep= "\t")

#### Join non-neutral table with KO function predictions
#non_neutral_ko_table <- inner_join(ko_function, sig_neutral_table, by ="id") 

#neutral_ko_name <- paste0(biosample,"_neutral_ko_contrib.tsv.gz")
#write.table(non_neutral_ko_table,gzfile(neutral_ko_name),sep = "\t")

print(paste("generating above path contributions across family and function"))
## Above pathway abundance
non_neutral_path_contrib %>% filter(model == "above" & padj <= 0.05) %>% group_by(function.) %>%
  summarize(genome_path_count = sum(genome_function_count), 
            mean_path_abun = mean(genome_function_count),
            genome_path_sd = sd(genome_function_count),
            unique_family_count = n_distinct(Phylum)) -> above_path_abund

above_path_name <- paste0(biosample,"_above_pathway_Phylum_table.tsv")
write.table(above_path_abund, file =above_path_name,sep="\t",row.names = FALSE, col.names=TRUE)

## Above family abundance
non_neutral_path_contrib %>% filter(model == "above" & padj <= 0.05) %>% group_by(Phylum) %>%
  summarize(taxon_path_count = sum(taxon_function_abun), 
            tax_path_mean_abun = mean(taxon_function_abun),
            tax_path_sd = sd(taxon_function_abun), 
            unique_function_count = n_distinct(function.)) -> above_family_abund

above_family_name <- paste0(biosample,"_above_pathway_Phylum_count.tsv")
write.table(above_family_abund, file =above_family_name,sep="\t",row.names = FALSE, col.names=TRUE)

print(paste("generating below path contributions"))
## Below pathway abundance
non_neutral_path_contrib %>% filter(model == "below" & padj <= 0.05) %>% group_by(function.) %>%
  summarize(genome_path_count = sum(genome_function_count), 
            mean_path_abun = mean(genome_function_count),
            genome_path_sd = sd(genome_function_count),
            unique_family_count = n_distinct(Phylum)) -> below_path_abund

below_pathway_name <- paste0(biosample,"_below_pathway_Phylum_table.tsv")
write.table(below_path_abund, file =below_pathway_name,sep="\t",row.names = FALSE, col.names=TRUE)

#below family
non_neutral_path_contrib %>% filter(model == "below" & padj <= 0.05) %>% group_by(Phylum) %>%
  summarize(taxon_path_count = sum(taxon_function_abun), 
            tax_path_mean_abun = mean(taxon_function_abun),
            tax_path_sd = sd(taxon_function_abun), 
            unique_function_count = n_distinct(function.)) -> below_family_abund

below_family_name <- paste0(biosample,"_below_pathway_Phylum_count.tsv")
write.table(below_family_abund, file =below_family_name,sep="\t",row.names = FALSE, col.names=TRUE)

print(paste("Finished"))