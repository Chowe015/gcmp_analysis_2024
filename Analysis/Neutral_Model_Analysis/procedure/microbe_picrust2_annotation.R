#load libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

## Step 5: Uses ggpicrust2 to run comparative analysis from combined psudo
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
comp <-args[1]

sink_name = paste0(comp,".","picrust2_annotation_log.txt")
sink(sink_name,append=FALSE,split=TRUE)
# Load metadata as a tibble
map_name <- paste0("./",comp,"/",comp,"_combined_metadata.txt")
metadata <- read_delim(map_name,delim = "\t", escape_double = FALSE,trim_ws = TRUE)
colnames(metadata)[1] <-"sample_name"

# Load KEGG pathway abundance
kegg_name <- paste0("./",comp,"_picrust2_out_pipeline","/",comp,"_KO_pred_metagenome_unstrat.tsv")
kegg_abundance <- ko2kegg_abundance(kegg_name)

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "neutrality", daa_method = "ALDEx2",reference = "neutral")

# Filter results for ALDEx2_Welch's t test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating KO to Kegg descriptions..."))
# Annotate pathway results using KO to KEGG conversion
kegg_annotated <-pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE, organism ="eco")

# Print out Results tables
kegg_result_name <- paste0(comp,"_","_ALDEx2_glm_KO_DE_descrip.tsv")
write.table(kegg_annotated, file =kegg_result_name,sep="\t", row.names = FALSE)


# Workflow for MetaCyc Pathway and EC

# Load MetaCyc pathway abundance and metadata
path_name <- paste0("./",comp,"_picrust2_out_pipeline","/",comp,"_path_abun_unstrat.tsv")
print(paste("Importing ",path_name))
metacyc_table <- read.table(path_name,header = T,check.names = FALSE)

# Perform pathway DAA using ALDEx2 method
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_table %>% column_to_rownames("pathway"), metadata = metadata, group = "neutrality", daa_method = "ALDEx2",reference = "neutral")

## Filter methods
metacyc_daa_results_df <- metacyc_daa_results_df[metacyc_daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating Pathway descriptions..."))
# Annotate the results
annotated_metacyc_daa_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = FALSE,organism ="eco")

# Print out Results tables
metacyc_result_name <- paste0(comp,"_","ALDEx2_glm_metacyc_DE_descrip.tsv")
write.table(annotated_metacyc_daa_results_df, file =metacyc_result_name,sep="\t", row.names = FALSE)

# Workflow for MetaCyc Pathway and EC
# Load MetaCyc pathway abundance and metadata
EC_path_name <- paste0("./",comp,"_picrust2_out_pipeline","/",comp,"_EC_pred_metagenome_unstrat.tsv")
print(paste("Importing",EC_path_name))
EC_table <- read.delim(EC_path_name,header = T,check.names = FALSE)

# Perform pathway DAA using LinDA method
EC_daa_results_df <- pathway_daa(abundance = EC_table %>% column_to_rownames("function"), metadata = metadata, group = "neutrality", daa_method = "ALDEx2",reference = "neutral")

## Filter methods
EC_daa_results_df <- EC_daa_results_df[EC_daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating EC Enzyme pathway Description"))
# Annotate the results
annotated_EC_daa_results_df <- pathway_annotation(
  pathway = "EC",
  daa_results_df = EC_daa_results_df,
  ko_to_kegg = FALSE,organism ="eco")

# Print out Results tables
EC_result_name <- paste0(comp,"_","ALDEx2_glm_EC_DE_descrip.tsv")
write.table(annotated_EC_daa_results_df, file =EC_result_name,sep="\t", row.names = FALSE)
print(paste("finished"))