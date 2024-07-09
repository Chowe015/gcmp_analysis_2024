#Loading Required libraries

library(qiime2R)
library(phyloseq)
library(tidyverse)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

feature_table_path <-args[1]
print(paste("feature_table",feature_table_path))

metadata_path <-args[2]
taxonomy_path <-args[3]

tree_path <-args[4]
print(paste("tree_path",tree_path))

#Import from .qza file into a phyloseq object
asv <- qza_to_phyloseq(features = feature_table_path)

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names=1, sep="\t")

### Import Tree file from biom output tree.nwk
tree <- read_tree(tree_path)

### Import taxonomy from biom output as .tsv format using read.table
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)


print(paste("metadata_path",metadata_path))
print(paste("taxonomy_path",taxonomy_path))
print(paste("Loading metadata files from path:", metadata_path))
print(paste("Loading feature table from path:", feature_table_path))
print(paste("Loading phylogenetic tree file from path:", tree_path))
print(paste("Loading Taxonomy text files from path:", taxonomy_path))

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
SAMPLE <- sample_data(metadata)

# Create Working phyloseq object
main_phylo <- phyloseq(OTU,tax1,SAMPLE,tree)

# visualize tax table
#table(tax_table(main_coral)[,"Kingdom"])

# remove unknown bacteria or unassigned
phylo_noking <-main_phylo %>%
  phyloseq::subset_taxa(!Kingdom %in% c("Unassigned","Unclassified d__Bacteria","d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Phylum %in% c("Unclassified Unassigned","Unclassified d__Archaea","Unclassified d__Bacteria","Unclassified d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Family %in% c("Mitochondria","Chloroplast"))

# verify present taxa in table
table(tax_table(main_phylo)[,"Kingdom"])

# Remove any samples with less than 1 read
phylo = prune_samples(sample_sums(phylo_noking)>1, phylo_noking)
phylo

# Subset only corals from database
subject <-subset_samples(phylo, sample_type_EMP == "coral" & outgroup == "n")# & within_group == "yes") # all corals without out groups, environmental samples & subset expeditions (i.e, singapore).
subject

print(paste("Generating Rarefied Coral dataset..."))
# rarefy to even depth
rarefied = rarefy_even_depth(subject, rngseed=111, sample.size=2000, replace=F, trimOTUs = TRUE)
rarefied

paste(print("Joining OTU and taxonomy tables from Rarefied Phyloseq object..."))
# Join otu and taxonomy tables by id
phyloseq::tax_table(rarefied)%>%
        as.data.frame()%>%
        rownames_to_column("id")%>%
        right_join(phyloseq::otu_table(rarefied)%>%
        as.data.frame()%>%
        rownames_to_column("id")) -> rare_otu_table

paste(print("Printing Rarefied ASV Table"))
## Output .csv from the biom file
write.csv(rare_otu_table, file ="phyloseq_rarefied_table_taxonomy.csv",row.names = FALSE)


print(paste("Agglomerate Taxonomy to the Family Level"))
##Agglomerate taxa to family 
glom <- tax_glom(rarefied, taxrank = 'Family', NArm = TRUE)

# Join otu and taxonomy tables by id
paste(print("Joining OTU and taxonomy tables from Agglomerated phyloseq object..."))
phyloseq::tax_table(glom)%>%
        as.data.frame() %>%
        rownames_to_column("id") -> glom_taxonomy

## Output .csv from the taxonomy file
write.csv(glom_taxonomy, file ="agglomerated_taxonomy.csv", row.names = FALSE)

print(paste("Generating Agglomerated ASV Table dataset..."))
phyloseq::otu_table(glom)%>%
        as.data.frame()%>%
        rownames_to_column("id") -> glom_otu_table

## Output .csv from the otu table file
write.csv(glom_otu_table, file ="agglomerated_table.csv",row.names = FALSE)

print(paste("Finished!"))
