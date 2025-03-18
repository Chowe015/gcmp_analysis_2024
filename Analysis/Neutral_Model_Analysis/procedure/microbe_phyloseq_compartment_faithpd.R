#Loading Required libraries

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(btools)
library(picante)
sink("Phyloseq_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
coral_tree_path <-args[5]
biosample <-args[6]
 

#Import from .qza file into a phyloseq object
asv <- qza_to_phyloseq(features = feature_table_path)
print(paste("feature_table",feature_table_path))

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names=1, sep="\t")

### Import Tree file from biom output tree.nwk
tree <- read_tree(tree_path)
print(paste("tree_path",tree_path))


### Import taxonomy from biom output as .tsv format using read.table
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

coral_tree <-read_tree(coral_tree_path)

print(paste("metadata_path",metadata_path))
print(paste("taxonomy_path",taxonomy_path))
print(paste("Loading feature table from path:", feature_table_path))
print(paste("Loading phylogenetic tree file from path:", tree_path))
print(paste("Loading Taxonomy text files from path:", taxonomy_path))
print(paste("Loading Coral Tree from path:", coral_tree_path))

##code referenced from Yan Hui: email me@yanh.org github: yanhui09
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}


# create matrix format for OTU and taxonomy table
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(tax.clean))


# Set metadata
print(paste("Loading metadata files from path:", metadata_path))
SAMPLE <- sample_data(metadata)

# Create Working phyloseq object
main_phylo <- phyloseq(OTU,tax1,SAMPLE,tree)


# visualize tax table
table(tax_table(main_phylo)[,"Kingdom"])

# remove unknown bacteria or unassigned
phylo_noking <-main_phylo %>%
  phyloseq::subset_taxa(!Kingdom %in% c("Unassigned","Unclassified d__Bacteria","d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Phylum %in% c("Unclassified Unassigned","Unclassified d__Archaea","Unclassified d__Bacteria","Unclassified d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Family %in% c("Mitochondria","Chloroplast"))

# verify present taxa in table
table(tax_table(phylo_noking)[,"Kingdom"])

# Remove any samples with less than 1 read
phylo = prune_samples(sample_sums(phylo_noking)>1, phylo_noking)
phylo


# Subset only corals from database
subject <-subset_samples(phylo, outgroup == "n" & tissue_compartment == print(paste(biosample)))
subject


# rarefy to even depth
print(paste("Generating Rarefied Coral dataset..."))
rarefied = rarefy_even_depth(subject, rngseed=111, sample.size=1000, replace=F, trimOTUs = TRUE)
rarefied

## Calculate Faiths Pd
print(paste("Calculating Faith Pd for host and microbiome"))

# pull out sample data 
phyloseq::sample_data(rarefied) %>%  
  group_by(sample_name_backup, expedition_number,BiologicalMatter) %>%
  as.data.frame() %>%
  select(sample_name_backup, expedition_number,BiologicalMatter) -> microbial_faith

#Microbial Phyloseq Pd analysis 
estimate_pd(rarefied) %>%
  as.data.frame() -> microbial_faith_pd

microbial_faith$faith_pd <- microbial_faith_pd$PD
microbial_faith$faith_SR <- microbial_faith_pd$SR


# Create data frame from sample data 
phyloseq::sample_data(rarefied) %>%  
  group_by(expedition_number, BiologicalMatter,Huang_Roy_tree_name) %>%
  as.data.frame() %>%
  select(expedition_number, BiologicalMatter,Huang_Roy_tree_name)-> test_df

# Create a new column titled eco to join expedition and biological matter
test_df$eco <-paste(test_df$expedition_number,test_df$BiologicalMatter, sep = "_")

#group and count total for each unique group maintaining NA values
test_df %>% group_by(eco, Huang_Roy_tree_name) %>%  summarise(counts=n()) %>%
  ungroup %>%
  complete(nesting(eco),
           nesting(Huang_Roy_tree_name),
           fill = list(quantity = 0)) -> test_table
# Fill NA with 0
test_table[is.na(test_table)] <-0

#Build Matrix
e <- unique(test_table$eco) 
t <- unique(test_table$Huang_Roy_tree_name) 
c <- test_table$counts

test_matrix <- matrix(c, nrow = length(e), ncol = length(t), byrow=TRUE)
rownames(test_matrix) = e
colnames(test_matrix) = t

# clean data set to match each other
clean_tree <- match.phylo.comm(phy = coral_tree, comm = test_matrix)$phy
clean_comm <- match.phylo.comm(phy = coral_tree, comm = test_matrix)$comm

coral_faith_pd <- pd(clean_comm, clean_tree, include.root=TRUE)
coral_faithpd_reorded <-coral_faith_pd[order(coral_faith_pd$PD, decreasing=TRUE),] 

write.table(microbial_faith, file =paste0(biosample,"_","microbial_faithpd_table.csv"), sep = ",",row.names =TRUE, col.names = TRUE)

write.csv(coral_faithpd_reorded, file =paste0(biosample,"_","Host_faithpd_table.csv") ,row.names = TRUE)

#### create ASV tables by id ** This file will be used in microbe_neutral_compartment.R and picrust2_neutral_table_generator.R

print(paste("Generating Agglomerated ASV Table dataset..."))
phyloseq::otu_table(glom)%>%
  as.data.frame()%>%
  rownames_to_column("id") -> glom_otu_table

## Output .csv from the otu table file
glom_otu_name <- paste0(biosample,"_glom_table.csv")
write.csv(glom_otu_table, file =glom_otu_name ,row.names = FALSE)

#### create taxonomy tables by id  ** This 

paste(print("Printing Agglomerated Taxonomy Table"))
phyloseq::tax_table(glom)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> glom_taxonomy

## Output .csv from the taxonomy file
glom_taxonomy_name <- paste0(biosample,"_glom_taxonomy.csv")
write.csv(glom_taxonomy, file =glom_taxonomy_name, row.names = FALSE)

#### Creating metadata file for downstream analysis ** This file is used as a mapping file for comparartive analysis for subset datasets MST. 

paste(print("Creating glom mapping file"))
phyloseq::sample_data(glom)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> glom_mapping

## Output .csv from the taxonomy file
glom_mapping_name <- paste0(biosample,"_glom_metadata.csv")
write.csv(glom_mapping, file =glom_mapping_name, row.names = FALSE)

#### Subset taxonomy tables ** This table contains non-agglomerated ASV ID which can be used 
#### for comparative analysis of significant non-neutral microbes. 

## Subset Rarefied taxonomy 
paste(print("Subset taxonomy rarefied phyloseq object..."))
phyloseq::tax_table(rarefied)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> rare_taxonomy

## Output .csv from the rarefied taxonomy file
rare_file_name <- paste0(biosample,"_rarefied_taxonomy.csv")
write.csv(rare_taxonomy, file =taxonomy_file_name, row.names = FALSE)

## Subset Rarefied otu table 
paste(print("Subset taxonomy rarefied phyloseq object..."))
phyloseq::otu_table(rarefied)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> rare_otu_table

## Output .csv from the rarefied taxonomy file
rare_otu_name <- paste0(biosample,"_rarefied_table.csv")
write.csv(rare_otu_table, file =rare_otu_name, row.names = FALSE)




print(paste("Finished!"))
