#Loading Required libraries

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(btools)
library(picante)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(minpack.lm)
library(Hmisc)
library(stats4)

sink("Phyloseq_neutral_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
#coral_tree_path <-args[5]

####Import from .qza file into a phyloseq object

print(paste("feature_table",feature_table_path))
asv <- qza_to_phyloseq(features = feature_table_path)

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names = 1, sep="\t")

#### Import Tree file from biom output tree.nwk

print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)

#coral_tree <-read_tree(coral_tree_path)

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

for (v in 1:nrow(tax.clean)){
  if (tax.clean[v,7] != ""){
    tax.clean$Species[v] <- paste(tax.clean$Genus[v], tax.clean$Species[v], sep = " ")
  } else if (tax.clean[v,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[v,1], sep = " ")
    tax.clean[v, 2:7] <- kingdom
  } else if (tax.clean[v,3] == ""){
    phylum <- paste("Unclassified", tax.clean[v,2], sep = " ")
    tax.clean[v, 3:7] <- phylum
  } else if (tax.clean[v,4] == ""){
    class <- paste("Unclassified", tax.clean[v,3], sep = " ")
    tax.clean[v, 4:7] <- class
  } else if (tax.clean[v,5] == ""){
    order <- paste("Unclassified", tax.clean[v,4], sep = " ")
    tax.clean[v, 5:7] <- order
  } else if (tax.clean[v,6] == ""){
    family <- paste("Unclassified", tax.clean[v,5], sep = " ")
    tax.clean[v, 6:7] <- family
  } else if (tax.clean[v,7] == ""){
    tax.clean$Species[v] <- paste("Unclassified ",tax.clean$Genus[v], sep = " ")
  }
}


### create matrix format for OTU and taxonomy table

print(paste("Loading metadata files from path:", metadata_path))
# create matrix format for OTU and taxonomy table
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(tax.clean))

# Set metadata
SAMPLE <- sample_data(metadata)
sample_names(SAMPLE) <- sample_names(OTU)

# Create Working phyloseq object
main_phylo <- phyloseq(OTU,tax1,SAMPLE,tree)
main_phylo

# visualize tax table
table(tax_table(main_phylo)[,"Kingdom"])

# Remove any samples with less than 1 read
phylo = prune_samples(sample_sums(main_phylo)>1, main_phylo)
phylo

# Calculate DISTANCE MATRICES
## Mucus

bray_dis = phyloseq::distance(phylo, method="bray")
uni_dis = phyloseq::distance(phylo, method="unifrac")
wuni_dis= phyloseq::distance(phylo, method="wUnifrac")

# Create data frame for downstream analysis
muc_df <- data.frame(sample_data(phylo))


# Ordination from Distance Matrix and PCoA Plots 

##Calculate the PCoA on Bray-Curtis Corals
rt.pcoa = ordinate(phylo, method="PCoA", distance=bray_dis)
# Set variables to zero for subsequent changes
#Ofra
pcoa<-0
pcoa <- plot_ordination(phylo, rt.pcoa ,color="neutraility")+  geom_point(size=3) + theme_bw() + ggtitle("Neutral Microbe Communities Bray-Curtis Distance")+theme(plot.title = element_text(hjust = 0.5)) + stat_ellipse(level = 0.95,type = "norm",aes(group=neutraility)) + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") + theme(text = element_text(size = 20,face = "bold"))) 

##Calculate the PCoA on Unifrac Corals
rt.pcoa2 = ordinate(phylo, method="PCoA", distance=uni_dis)
# Set variables to zero for subsequent changes
#Ofra
pcoa2<-0
pcoa2 <- plot_ordination(phylo, rt.pcoa2 ,color="neutraility")+  geom_point(size=3) + theme_bw() + ggtitle("Neutral Microbe Communities UNifrac Distance")+theme(plot.title = element_text(hjust = 0.5)) + stat_ellipse(level = 0.95,type = "norm",aes(group=neutraility)) + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") + theme(text = element_text(size = 20,face = "bold"))) 

##Calculate the PCoA on Unifrac Corals
rt.pcoa3 = ordinate(phylo, method="PCoA", distance=wuni_dis)
# Set variables to zero for subsequent changes
#Ofra
pcoa3<-0
pcoa3 <- plot_ordination(phylo, rt.pcoa3 ,color="neutraility")+  geom_point(size=3) + theme_bw() + ggtitle("Neutral Microbe Communities \n Weighted UNifrac Distance")+theme(plot.title = element_text(hjust = 0.5)) + stat_ellipse(level = 0.95,type = "norm",aes(group=neutraility)) + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") + theme(text = element_text(size = 20,face = "bold"))) 

print(paste("Saving Neutral PcoA Graphs!"))
neutral_bray_name <- paste0(biosample,"_neutrality_bray-curtis.pdf")
ggsave(neutral_dot, filename=neutral_dotplot_name, width = 20, height = 10, dpi= 600)



print("finished!")