#Loading Required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(minpack.lm)
library(Hmisc)
library(stats4)

sink("merged_tables_log.txt",append=FALSE,split=TRUE)


#Get user input and assign to variables
setwd("C:/Users/Colin Howe/Desktop/PSU Doctoral Program/GCMP_Working")

## Merging psudo_sample tables taxonomy and mapping files across compartments
#Mucus
mucus_table <-read.csv("./Mucus/Mucus_combinded_table.csv", row.names=1, check.names=FALSE)
mucus_tax <-read.csv("./Mucus/Mucus_qiime_neutral_taxonomy.csv", row.names=1, check.names=FALSE)
mucus_mapping <-read.csv("./Mucus/Mucus_qiime_neutral_mapping.csv", row.names=1, check.names=FALSE)

#Tissue
tissue_table <-read.csv("./Tissue/Tissue_combinded_table.csv", row.names=1, check.names=FALSE)
tissue_tax <-read.csv("./Tissue/Tissue_qiime_neutral_taxonomy.csv", row.names=1, check.names=FALSE)
tissue_mapping <-read.csv("./Tissue/Tissue_qiime_neutral_mapping.csv", row.names=1, check.names=FALSE)

#Skeleton
skeleton_table <-read.csv("./Skeleton/Skeleton_combinded_table.csv", row.names=1, check.names=FALSE)
skeleton_tax <-read.csv("./Skeleton/Skeleton_qiime_neutral_taxonomy.csv", row.names=1, check.names=FALSE)
skeleton_mapping <-read.csv("./Skeleton/Skeleton_qiime_neutral_mapping.csv", row.names=1, check.names=FALSE)


## Merge OTU Tables across compartments
merge(mucus_table,tissue_table, by=0, all=T)%>% 
  merge(skeleton_table) -> merge_table 

## Assign all NA as 0
merge_table[is.na(merge_table)] <-0

## Rename the first column name to blank
colnames(merge_table)[1] <- ""

## Print output files
print(paste("Writing Combined Psudo_table across compartments"))
psudo_table_name <- paste0("combined_psudo_table.tsv.gz")
write.table(merge_table, gzfile(psudo_table_name), sep= "\t")

## Merge Taxonomy files 
merge(mucus_tax,tissue_tax, by=0, all=T) %>% 
  merge(skeleton_tax) -> merge_tax 

## coalesce taxonomy columnss
merge_tax %>% mutate(
  Taxon = coalesce(Taxon.x,Taxon.y),
  Consensus = coalesce(Consensus.x,Consensus.y)) %>%
  select(Row.names,Taxon,Consensus) -> combined_tax

## Rename the First column name to blank
colnames(combined_tax)[1] <- ""

# Print output files
print(paste("Writing Combined Psudo_taxonomy across compartments"))
psudo_tax_name <- paste0("combined_psudo_tax.tsv.gz")
write.table(combined_tax,gzfile(psudo_tax_name), sep= "\t")

## Merging Mapping Files
## Filter row artifacts
m_map <-mucus_mapping[-c(1,2,3),]
t_map <-tissue_mapping[-c(1,2,3),]
s_map <-skeleton_mapping[-c(1,2,3),]

## Create new column with compartment and neutrality
s_map %>% mutate(compartment_neutrality  = paste(compartment, model, sep = '_')) -> skeleton_map
t_map %>% mutate(compartment_neutrality  = paste(compartment, model, sep = '_')) -> tissue_map
m_map %>% mutate(compartment_neutrality  = paste(compartment, model, sep = '_')) -> mucus_map

## Merge each mapping file
merge(mucus_map,tissue_map, by=0, all=T) %>% 
  merge(skeleton_map) -> merge_map 

## Coalesce columns
merge_map %>% mutate(compartment = coalesce(compartment.x,compartment.y),
                     model = coalesce(model.x,model.y),
                     compartment_neutrality = coalesce(compartment_neutrality.x, compartment_neutrality.y)) %>% 
  select(Row.names,model,compartment,compartment_neutrality)-> combined_map

#rename column Row.names
colnames(combined_map)[1] <- ""

# Print output files
print(paste("Writing Combined Psudo_mapping across compartments"))
psudo_map_name <- paste0("combined_psudo_mapping.tsv")
write.table(combined_map,psudo_map_name,row.names = FALSE, col.names = TRUE)




print(paste("Finished!"))
