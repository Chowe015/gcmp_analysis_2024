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

comp <- c("s_disp_m","s_disp_t","t_disp_m")

for (i in comp){
sink_name = paste0(i,"_","picrust2_annotation_log.txt")
sink(sink_name,append=FALSE,split=TRUE)

# Load metadata as a tibble
map_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/",i,"_combined_metadata.txt")
metadata <- read_delim(map_name,delim = "\t", escape_double = FALSE,trim_ws = TRUE)
colnames(metadata)[1] <-"sample_name"

metadata %>% filter(neutrality=="neutral") -> map

print(paste0("***** processing ", i," dispersal ******"))

# Load KEGG pathway abundance
kegg_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv")
kegg_abundance <- ko2kegg_abundance(kegg_name)

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = map, group = "neutrality", daa_method = "ALDEx2",reference = "below")

# Filter results for ALDEx2_Welch's t test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating KO to Kegg descriptions..."))
# Annotate pathway results using KO to KEGG conversion
kegg_annotated <-pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

# Print out Results tables
kegg_result_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/KO_metagenome_out/",i,"_ALDEx2_glm_KO_DE_descrip.tsv")
write.table(kegg_annotated, file =kegg_result_name,sep="\t", row.names = FALSE)


# Workflow for MetaCyc Pathway and EC
# Load MetaCyc pathway abundance and metadata
path_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv")
print(paste("Importing ",path_name))
metacyc_table <- read.table(path_name,header = T,check.names = FALSE)

# Perform pathway DAA using ALDEx2 method
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_table %>% column_to_rownames("pathway"), metadata = map, group = "neutrality", daa_method = "ALDEx2",reference = "below")

## Filter methods
metacyc_daa_results_df <- metacyc_daa_results_df[metacyc_daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating Pathway descriptions..."))
# Annotate the results
annotated_metacyc_daa_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = FALSE)

# Print out Results tables
metacyc_result_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/pathway_out/",i,"_ALDEx2_glm_metacyc_DE_descrip.tsv")
write.table(annotated_metacyc_daa_results_df, file =metacyc_result_name,sep="\t", row.names = FALSE)

# Workflow for MetaCyc Pathway and EC
# Load MetaCyc pathway abundance and metadata
EC_path_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv")
print(paste("Importing",EC_path_name))
EC_table <- read.delim(EC_path_name,header = T,check.names = FALSE)

# Perform pathway DAA using LinDA method
EC_daa_results_df <- pathway_daa(abundance = EC_table %>% column_to_rownames("function"), metadata = map, group = "neutrality", daa_method = "ALDEx2",reference = "below")

## Filter methods
EC_daa_results_df <- EC_daa_results_df[EC_daa_results_df$method == "ALDEx2_glm test", ]

print(paste("Annotating EC Enzyme pathway Description"))
# Annotate the results
annotated_EC_daa_results_df <- pathway_annotation(
  pathway = "EC",
  daa_results_df = EC_daa_results_df,
  ko_to_kegg = FALSE)

# Print out Results tables
EC_result_name <- paste0("./final_picrust2/",i,"_picrust2_out_pipeline/EC_metagenome_out/",i,"_ALDEx2_glm_EC_DE_descrip.tsv")
write.table(annotated_EC_daa_results_df, file =EC_result_name,sep="\t", row.names = FALSE)

}

print(paste("finished"))