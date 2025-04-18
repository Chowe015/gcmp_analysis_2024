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

sink("Picrust2_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

print(paste("importing data..."))
feature_table_path <-args[1]
neutral_results_path <-args[2]
kegg_path <-args[3]
pathway_path <-args[4]
biosample <-args[5]

#### Read glom and taxonomy input data
asv_table <- read.csv(feature_table_path, header =T, sep = ",", check.names=FALSE) 

#### Read neutral table input
neutral_results <- read.csv(neutral_results_path, header=T,comment.char="", sep=",", check.names = FALSE)

#### Join neutral and glom tables include non-significant taxa
print(paste("Joining neutral results & Glom  table"))
inner_join(neutral_results,asv_table, by ="id") %>% select(-c(freq.pred,pred.lwr,pred.upr,bino.pred,bino.lwr,bino.upr)) -> neutral_tax_table

## Write neutral tax table
neutral_tax_table_name <- paste0(biosample,"all_neutral_tax_table.csv")
write.table(neutral_tax_table, file =neutral_tax_table,sep=",",row.names = FALSE, col.names=TRUE)

## Filter for significant families
print(paste("Filtering Significant ASV from Neutral Table"))
filter(neutral_tax_table, padj < 0.05) %>% select(c(id,padj,Phylum,Class,Order,Family,model)) -> sig_neutral_table

#### Inner_join feature table with significant microbial families 
print(paste("joining neutral table and ref taxonomy"))
nonneutral_tax <- inner_join(asv_table,sig_neutral_table, by = c("id","Phylum","Class","Order","Family"))

sig_tax_table_name <- paste0(biosample,"_sig_neutral_tax_table.csv")
write.table(nonneutral_tax, file =sig_tax_table_name,sep=",",row.names = FALSE, col.names=TRUE)

## Maybe I need to merge the neutral model results with the rarefied dataset to have asv level tables of only the sig_non_neutral_families

print(paste("Importing KO and Kegg contribution tables"))
pathway <- read.table(gzfile(pathway_path), header =T, sep = "\t")#,check.names = FALSE)
colnames(pathway)[which(names(pathway)== "taxon")] <- "id"

ko_function <- read.table(gzfile(kegg_path), header =T, sep = "\t")#, check.names = FALSE)
colnames(ko_function)[which(names(pathway)== "taxon")] <- "id"

print(paste("joining neutral table with pathway & KO data tables"))
#### Join non-neutral table with pathway predictions
non_neutral_pathway_contrib <- inner_join(pathway,nonneutral_tax, by ="id") 

#### Join non-neutral table with KO function predictions
non_neutral_ko_table <- inner_join(ko_function, nonneutral_tax, by ="id") 

print(paste("Output Merged Tables"))

# Output Tables
neutral_table_name <- paste0(biosample,"_sig_nonneutral_feature_table.csv")
write.table(nonneutral_tax, file =neutral_table_name,sep=",",row.names = FALSE, col.names=TRUE)

neutral_pathway_name <- paste0(biosample,"_neutral_pathway_contrib.tsv.gz")
write.table(non_neutral_pathway_contrib, gzfile(neutral_pathway_name),sep= "\t")

neutral_ko_name <- paste0(biosample,"_neutral_ko_contrib.tsv.gz")
write.table(non_neutral_ko_table,gzfile(neutral_ko_name),sep = "\t")

## Above pathway abundance
non_neutral_pathway_contrib %>% filter(model == "above") %>% group_by(function.) %>%
  summarize(genome_path_count = sum(genome_function_count), 
            mean_path_abun = mean(genome_function_count),
            genome_path_sd = sd(genome_function_count),
            unique_family_count = n_distinct(Family)) -> above_path_abund

above_path_name <- paste0(biosample,"_above_pathway_table.csv")
write.table(above_path_abund, file =above_path_name,sep=",",row.names = FALSE, col.names=TRUE)

## Above family abundance
non_neutral_pathway_contrib %>% filter(model == "above") %>% group_by(Family) %>%
  summarize(taxon_path_count = sum(taxon_function_abun), 
            tax_path_mean_abun = mean(taxon_function_abun),
            tax_path_sd = sd(taxon_function_abun), 
            unique_function_count = n_distinct(function.)) -> above_family_abund

above_family_name <- paste0(biosample,"_above_family_table.csv")
write.table(above_family_abund, file =above_family_name,sep=",",row.names = FALSE, col.names=TRUE)

## Below pathway abundance
non_neutral_pathway_contrib %>% filter(model == "below") %>% group_by(function.) %>%
  summarize(genome_path_count = sum(genome_function_count), 
            mean_path_abun = mean(genome_function_count),
            genome_path_sd = sd(genome_function_count),
            unique_family_count = n_distinct(Family)) -> below_path_abund

below_pathway_name <- paste0(biosample,"_below_pathway_table.csv")
write.table(below_path_abund, file =below_pathway_name,sep=",",row.names = FALSE, col.names=TRUE)

#below family
non_neutral_pathway_contrib %>% filter(model == "below") %>% group_by(Family) %>%
  summarize(taxon_path_count = sum(taxon_function_abun), 
            tax_path_mean_abun = mean(taxon_function_abun),
            tax_path_sd = sd(taxon_function_abun), 
            unique_function_count = n_distinct(function.)) -> below_family_abund

below_family_name <- paste0(biosample,"_below_family_table.csv")
write.table(below_family_abund, file =below_family_name,sep=",",row.names = FALSE, col.names=TRUE)

print(paste("Finished"))