#Loading Required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(minpack.lm)
library(Hmisc)
library(bbmle)
library(biomformat)


## Step 2 Run Sloan Neutral Model to produce neutral table and figures
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

neutral_results_path <-args[1]
sink_table_path <- args[2]
sink_map_path <- args[3]
source_table_path <- args[4]
source_map_path <- args[5]
sink_name <- args[6]
source_name <- args[7]

## merge T_glom_table and M_glom_table and filter gcmp_glom_metadata for only Mucus and tissue and then pair with neutral results table
sink_title = paste0(sink_name,"_",source_name,"neutral_model_results_log.txt")
sink(sink_title,append=FALSE,split=TRUE)

## Import compartment glom files
#m_table <- read.table("./mucus/m_glom_table.tsv", sep = "\t",header = TRUE,check.name=FALSE)
#s_table <- read.table("./skeleton/s_glom_table.tsv", sep = "\t",header = TRUE,check.name=FALSE)
## Import Mapping files
#m_map <- read.table("./mucus/m_glom_metadata.tsv", sep = "\t",header = TRUE,check.name=FALSE)
#s_map <- read.table("./skeleton/s_glom_metadata.tsv", sep = "\t",header = TRUE,check.name=FALSE)
## Import neutral model results
#neutral_results_padj <-read.table("./dispersal/mucus/M_sub_S_abund_neutral_table.tsv", sep = "\t",header = TRUE,check.name=FALSE)
#disp_name <-"s_disp_m"


print(paste("Imported GCMP OTU Table"))
## Import glom compartment table
sink_table <-read.table(sink_table_path, sep = "\t",header = TRUE,check.name=FALSE)
source_table <-read.table(source_table_path, sep = "\t",header = TRUE,check.name=FALSE)

##Import mapping files
sink_map <- read.table(sink_map_path, sep = "\t",header = TRUE,check.name=FALSE)
source_map <- read.table(source_map_path, sep = "\t",header = TRUE,check.name=FALSE)

## Import Neutral Model Results 
neutral_results_padj <-read.table(neutral_results_path, sep = "\t",header = TRUE,check.name=FALSE)

## Merge compartment tables together
full_join(sink_table, source_table, by ="id") -> merged_table

## Merge compartment tables together
rbind(sink_map, source_map) -> glom_mapping

############################### create compartment pseudo table outputs #########################
print(paste("Writing Above and Below Pseudo tables"))

# Subset Taxonomy and join with neutral model results
b <-"below"
neutral = c("above","below")
for (b in neutral){
    ## Inner_join will connect  ASV id with taxonomy using id column
  #inner_join(neutral_results_padj,glom_tax, by = "id") -> non_neutral_tax
  
  ## Subset ASV from each model and significance 
  sig_neutral_tax = select(filter(neutral_results_padj, model == b & padj <=0.05),id,p_abundance,freq,padj,Phylum,Class,Order,Family,model)
  
  #subset otu table based on neutral results and rename column names
  neutral_table <-subset(merged_table, id  %in%  sig_neutral_tax$id)
  l <- length(neutral_table)
  colnames(neutral_table)[2:l] <- paste(colnames(neutral_table)[2:l], b, sep = "_")
 colnames(neutral_table)[1] <- "id" 
  ## Create new columns for the metadata data
 glom_mapping %>% mutate(psuedo_sample_name = paste(colnames(neutral_table)[2:l]),
                          compartment_neutrality = paste(sink_name,source_name,b, sep = '_'),
                          compartment = paste (sink_name,source_name, sep = '_'),
                          neutrality = paste(b)) %>% 
    select(psuedo_sample_name,sample_name_backup, host_scientific_name,expedition_number,political_area,ocean_area,compartment_neutrality,compartment,neutrality) %>% 
    as.data.frame() -> glom_metadata
  
  ## Subset significant results with taxonomy to retain only significant taxonomy
  #subset_tax <- subset(glom_tax, id %in% sig_neutral_tax$id) 
  #psuedo_tax_name <- paste0(disp_name,"_",b,"_taxonomy.tsv")
  #(write.table(subset_tax, file =psuedo_tax_name,sep="\t", row.names = FALSE))
  
  psuedo_table_name <- paste0("./pseudo/",sink_name,"_",source_name,"_",b,"_table.tsv")
  (write.table(neutral_table, file =psuedo_table_name,sep="\t", row.names = FALSE))
  
  psuedo_map_name <- paste0("./pseudo/",sink_name,"_",source_name,"_",b,"_metadata.txt")
  (write.table(glom_metadata, file =psuedo_map_name,sep="\t", row.names = FALSE))
  
}

########### Redo for neutral without filtering for  significant padj ########
print(paste("Writing non-significant taxa table"))
## Subset ASV from each model and signifcance 
sig_neutral_ref = select(filter(neutral_results_padj, model == "neutral" & padj >=0.05),id,p_abundance,freq,padj,Phylum,Class,Order,Family,model)

#subset otu table based on neutral results and rename column names
neutral_table_test <-subset(merged_table, id  %in%  sig_neutral_ref$id)
la = length(neutral_table_test)
colnames(neutral_table_test)[2:la] <- paste(colnames(neutral_table_test)[2:la], "neutral", sep = "_")
colnames(neutral_table_test)[1] <- "id"

## Create new columns for the metadata data
glom_mapping %>% mutate(psuedo_sample_name = paste(colnames(neutral_table_test)[2:la]),
                            compartment_neutrality = paste(sink_name,source_name,"neutral", sep = '_'),
                            compartment = paste (sink_name,source_name, sep = '_'),
                            neutrality = paste("neutral")) %>% 
  select(psuedo_sample_name,sample_name_backup, host_scientific_name,expedition_number,political_area,ocean_area,compartment_neutrality,compartment,neutrality) %>% 
  as.data.frame() -> neutral_glom_metadata

## Subset significant results with taxonomy to retain only significant taxonomy
#subset_tax <- subset(glom_tax, id %in% sig_neutral_tax$id) 

## Output the tables for downstream analysis 
psuedo_table_name <- paste0("./pseudo/",sink_name,"_",source_name,"_","neutral_table.tsv")
(write.table(neutral_table_test, file =psuedo_table_name, sep="\t",row.names = FALSE))

psuedo_map_name <- paste0("./pseudo/",sink_name,"_",source_name,"_","neutral_metadata.txt")  ## Remember to combined mapping files in Excel
(write.table(neutral_glom_metadata, file =psuedo_map_name,sep="\t", row.names = FALSE))

################### Merge OTU Tables, Taxonomy & Metadata across neutrality ####################################

############### Import Tables ###############
above_import = paste0("./pseudo/",sink_name,"_",source_name,"_above_table.tsv")
below_import = paste0("./pseudo/",sink_name,"_",source_name,"_below_table.tsv")
neutral_import = paste0("./pseudo/",sink_name,"_",source_name,"_neutral_table.tsv")

above_table <- read.table(above_import,header = TRUE,sep="\t", check.names=FALSE)
neutral_table <- read.table(neutral_import,header = TRUE,sep="\t", check.names=FALSE)
below_table <- read.table(below_import,header = TRUE,sep="\t", check.names=FALSE)

## 1) Merge Tables first 
full_join(above_table, below_table, by ="id") -> merge_table

#merge_table %>% mutate(id = coalesce(id.x,id.y)) %>% 
#  relocate(id) %>%  
#  select(!c(id.x,id.y,Row.names))-> merge_test
#merge_table[is.na(merge_table)] <-0

full_join(merge_table,neutral_table, by ="id") -> full_table
#full_table %>% mutate(id = coalesce(id.x,id.y)) %>% 
#  relocate(id) %>%  
#  select(!c(id.x,id.y,Row.names))-> full_table
full_table[is.na(full_table)] <-0

colnames(full_table)[1] <- "#OTU ID"

## Print output files
print(paste("Writing Combined pseudo_table across compartments"))
pseudo_table_name <- paste0("./pseudo/",sink_name,"_",source_name,"_combined_pseudo_table.tsv")
write.table(full_table, file=pseudo_table_name, sep="\t",row.names = FALSE)

############### Merge Taxonomy files ###############
tax2_import = paste0("./taxonomy.tsv")
taxonomy2 <- read.table(tax2_import, sep = "\t",header = TRUE,check.name=FALSE)
colnames(taxonomy2)[1] <-"id"
neutral_tax <-subset(taxonomy2, id  %in%  full_table$id)
colnames(neutral_tax)[1] <-"Feature ID"

# Print output files
print(paste("Writing pseudo_taxonomy across compartments"))
pseudo_tax_name <- paste0("./pseudo/",sink_name,"_",source_name,"_pseudo_tax.tsv")
write.table(neutral_tax,pseudo_tax_name, sep= "\t", row.names = FALSE)


############### Merge Metadata Files ###############
print(paste("Joining",sink_name,"_",source_name,"mapping files"))

above_map_import = paste0("./pseudo/",sink_name,"_",source_name,"_above_metadata.txt")
below_map_import = paste0("./pseudo/",sink_name,"_",source_name,"_below_metadata.txt")
neutral_map_import = paste0("./pseudo/",sink_name,"_",source_name,"_neutral_metadata.txt")

above_map <- read.table(above_map_import,header = TRUE,sep="\t", check.names=FALSE)
below_map <- read.table(below_map_import,header = TRUE,sep="\t", check.names=FALSE)
neutral_map <- read.table(neutral_map_import,header = TRUE,sep="\t", check.names=FALSE)

rbind(above_map,below_map,neutral_map) -> combined_map

# Print output files
print(paste("Writing Combined Metadata File"))
comb_map_name <- paste0("./pseudo/",sink_name,"_",source_name,"_combined_metadata.txt")
write.table(combined_map,comb_map_name, sep= "\t", row.names = FALSE)

print(paste("Finished!"))

