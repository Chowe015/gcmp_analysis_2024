#Loading Required libraries

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(btools)

sink("Phyloseq_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
biosample <-args[5]

####Import from .qza file into a phyloseq object

print(paste("feature_table",feature_table_path))
asv <- qza_to_phyloseq(features = feature_table_path)

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names=1, sep="\t")

#### Import Tree file from biom output tree.nwk

print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)

#### Import taxonomy from biom output as .tsv format using read.table

print(paste("Loading Taxonomy text files from path:", taxonomy_path))
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

#        **code referenced from Yan Hui: email me@yanh.org github: yanhui09**

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


### create matrix format for OTU and taxonomy table

print(paste("Loading metadata files from path:", metadata_path))
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(tax.clean))


# Set metadata
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


### rarefy to even depth

print(paste("Generating Rarefied Coral dataset..."))
rarefied = rarefy_even_depth(subject, rngseed=111, sample.size=1000, replace=F, trimOTUs = TRUE)
rarefied

#### Agglomerate taxa to family 

print(paste("Agglomerated Taxonomy to the Family Level"))
glom <- tax_glom(rarefied, taxrank = 'Family', NArm = TRUE)

## Output .csv from the joined table

paste(print("Printing Glom ASV & taxononmy Table"))
rare_file_name <- paste0(biosample,"_","glom_table_taxonomy.csv")
write.csv(glom_otu_table, file =rare_file_name, row.names = FALSE)

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

#### Creating metadata file for downstream analysis ** This file is used as a mapping file for comparartive analysis for subset datasets MST. 

paste(print("Creating new subset mapping file verify samples are all from correct compartment"))
phyloseq::sample_data(glom)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> sampleID_table

## Output .csv from the taxonomy file
sample_file_name <- paste0(biosample,"_glom_metadata.csv")
write.csv(sampleID_table, file =sample_file_name, row.names = FALSE)


#### Creating metadata file for downstream analysis

paste(print("Creating glom mapping file"))
phyloseq::sample_data(glom)%>%
  as.data.frame() %>%
  rownames_to_column("id") -> glom_mapping

## Output .csv from the taxonomy file
glom_mapping_name <- paste0(biosample,"_glom_metadata.csv")
write.csv(glom_mapping, file =glom_mapping_name, row.names = FALSE)

##paste(print("Joining OTU and taxonomy tables from agglomerated Phyloseq object..."))
#### Join asv and taxonomy tables by id **This will be used in microbe_picrust2_neutral_table_generator.R
#phyloseq::tax_table(glom)%>%
#  as.data.frame()%>%
#  rownames_to_column("id")%>%
#  right_join(phyloseq::otu_table(glom)%>%
#               as.data.frame()%>%
#               rownames_to_column("id")) -> glom_otu_table
print(paste("Finished!"))
