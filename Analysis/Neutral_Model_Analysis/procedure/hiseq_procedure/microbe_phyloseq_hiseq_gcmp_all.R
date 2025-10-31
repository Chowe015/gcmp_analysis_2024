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
library(ggsignif)

#Step 1: Import qiime2 tables, mapping, taxonomy, tree and coral phylogeny across mucus tissue and skeleton
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
coral_tree_path <-args[5]

sink_name = paste0("step1_Phyloseq_Hiseq_Faith_PD_log.txt")
sink(sink_name,append=FALSE,split=TRUE)

many_col <- c("#A6CEE3", "black","#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00","#1B9E77","#B15928",
              "#FB9A99", "#CAB2D6", "#377EB8","#D95F02", "#7570B3","#FFFF99","#666666" ,"#66A61E", "#E6AB02","#A6761D", "#CAB2D6",
              "#B15928","#FFFF99", "#D95F02", "#E7298A", "#7570B3","#666666" , "#E6AB02", "brown3","#A6761D", "#377EB8", "#4DAF4A" ,
              "#984EA3", "#FFFF33","#FF7F00","cyan3", "#4DAF4A" ,"#FF7F00","#984EA3", "#A65628","#1B9E77", "#F781BF", "#F46D43","#ABDDA4",
              "#FDAE61","black", "#FEE08B","#E6F598","#3288BD","deeppink","#66C2A5","cyan1","red3", "#E69F00","#56B4E9","#999999","#F0E442",
              "#0072B2","#D55E00","#D53E4F","yellow","skyblue","coral4","#A6CEE3","white","darkkhaki","brown1","#FFFF33","chocolate",
              "darkorchid2","#FF7F00","#1B9E77","#66A61E","#FDBF6F", "#FB9A99","#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A",
              "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00","#1B9E77","#B15928", "#FB9A99", "#CAB2D6", "#377EB8","#D95F02", "#7570B3","#FFFF99",
              "#666666" ,"#66A61E", "#E6AB02")


####Import from .qza file into a phyloseq object

print(paste("feature_table",feature_table_path))
asv <- qza_to_phyloseq(features = feature_table_path)

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T,comment.char="", row.names=1, sep="\t")

#### Import Tree file from biom output tree.nwk

#print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)

coral_tree <-read_tree(coral_tree_path)

#### Import taxonomy from biom output as .tsv format using read.table

print(paste("Loading Taxonomy text files from path:", taxonomy_path))
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

#        **code referenced from Yan Hui: email me@yanh.org github: yanhui09**

tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
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
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(tax.clean))


# Set metadata
SAMPLE <- sample_data(metadata)

# create matrix format for OTU and taxonomy table
OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
tax1 = tax_table(as.matrix(tax.clean))

# Set metadata
SAMPLE <- sample_data(metadata)

# Create Working phyloseq object
main_phylo <- phyloseq(OTU,tax1,SAMPLE,tree)
main_phylo

# visualize tax table
table(tax_table(main_phylo)[,"Kingdom"])

# remove unknown bacteria or unassigned
phylo_noking <-main_phylo %>%
  phyloseq::subset_taxa(!Kingdom %in% c("Unassigned","Unclassified d__Bacteria","d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Phylum %in% c("Unclassified Unassigned","Unclassified d__Archaea","Unclassified d__Bacteria","Unclassified d__Eukaryota")) %>%
  phyloseq::subset_taxa(!Family %in% c("Mitochondria","Chloroplast"))

# visualize tax table
table(tax_table(phylo_noking)[,"Kingdom"])

# Subset only corals from database
subject <-subset_samples(main_phylo, outgroup == "n")
subject
print(paste("Agglomerating taxa to the Family level"))
###Agglomerate taxa to family 
glom <- tax_glom(subject, taxrank = 'Genus', NArm = TRUE)
glom

phyloseq::tax_table(glom)%>%
  as.data.frame()%>%
  rownames_to_column("id") -> glom_tax

print(paste("Calculating Alpha Diversity Results"))
gcmp_rich = estimate_richness(glom, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
comp.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, sample_data(glom)$BiologicalMatter, p.adjust.method = "fdr")
comp.shannon

region.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, sample_data(glom)$expedition_number, p.adjust.method = "fdr")
region.shannon

family.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, sample_data(glom)$family, p.adjust.method = "fdr")
family.shannon

## Treaments
alpha_compartment_all <- plot_richness(glom, x="BiologicalMatter", color= "BiologicalMatter", measures = "Shannon") +
  geom_boxplot(width=0.5) + geom_jitter()+ scale_color_manual(values=c("#009197","#8E3B97","#EC8C5C"))+ geom_point(na.rm = TRUE) +
  theme_classic() + theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") + geom_signif(comparisons = list(c("Coral Mucus", "Coral Tissue"), c("Coral Skeleton","Coral Mucus"),
                                                                                                                                 c("Coral Tissue", "Coral Skeleton")), y_position = c(6,5,5.5))

alpha_family_all <- plot_richness(glom, x="family", color= "family", measures = "Shannon") +
  geom_boxplot(width=0.5) + geom_jitter()+ scale_color_manual(values=many_col)+ geom_point(na.rm = TRUE) + theme_classic()+
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 90, vjust = 0, hjust = 0.5)) +
 theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") 

alpha_region_all <- plot_richness(glom, x="expedition_number", color= "expedition_number", measures = "Shannon") +
  geom_boxplot(width=0.5) + geom_jitter()+ scale_color_manual(values=many_col)+ geom_point(na.rm = TRUE) + theme_classic()+
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") 

#print(paste("Plotting Alpha Diversity Results"))
alpha_compartment_name <- paste0("Hiseq All Compartment Alpha Diversity.pdf")
ggsave(alpha_compartment_all, filename=alpha_compartment_name)

alpha_region_name <- paste0("Hiseq All GCMP Expedition Alpha Diversity.pdf")
ggsave(alpha_region_all, filename=alpha_region_name)

alpha_family_name <- paste0("Hiseq All Host Family Alpha Diversity.pdf")
ggsave(alpha_family_all, filename=alpha_family_name)


print(paste("Calculating bray_curtis Distance Matrices "))
## All samples
bray_curtis = phyloseq::distance(glom, method="bray")
print(paste("Calculating Unifrac Distance Matrices "))
unifrac = phyloseq::distance(glom, method="unifrac")
print(paste("Calculating weight_unifrac Distance Matrices "))
weight_unifrac = phyloseq::distance(glom, method="wUnifrac")

# Create data frame for downstream analysis
all_df <- data.frame(sample_data(glom))

print(paste("Calculating Beta Diversity Results"))
distance_methods <-c("bray_curtis","unifrac","weight_unifrac")

set.seed(129)
for (i in distance_methods){ 
  form <- as.formula(paste(i,"BiologicalMatter", sep="~"))
  print(form)
  adonis2(form, data=all_df)-> result
  print(result)
}

set.seed(129)
for (i in distance_methods){ 
  form <- as.formula(paste(i,"expedition_number", sep="~"))
  print(form)
  adonis2(form, data=all_df)-> result
  print(result)
}

set.seed(129)
for (i in distance_methods){ 
  form <- as.formula(paste(i,"family", sep="~"))
  print(form)
  adonis2(form, data=all_df)-> result
  print(result)
}

print(paste("Plotting Beta Diversity Results"))
##Calculate the PCoA on Compartment Region Family Corals

for (r in distance_methods){
print(paste0(r,"Ordinating Distance Matrix"))
rt.pcoa = ordinate(glom, method="PCoA", distance=r)

pcoa_title = paste0(r,"Distance Matrix for GCMP Hiseq 2018 Compartments")
pcoa <- plot_ordination(glom, rt.pcoa,color="BiologicalMatter") + scale_color_manual(values=c("#009197","#8E3B97","#EC8C5C"))+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=BiologicalMatter)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title)+theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 10))

pcoa1_title = paste0(r,"Distance Matrix for GCMP Hiseq 2018 Expeditions")
pcoa1 <- plot_ordination(glom, rt.pcoa,color="expedition_number") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=expedition_number)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa1_title)+theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 10))

pcoa2_title = paste0(r,"Distance Matrix for GCMP Hiseq 2018 Host Family")
pcoa2 <- plot_ordination(glom, rt.pcoa,color="family") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=family)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa2_title)+theme(plot.title = element_text(hjust = 0.5))  + theme(plot.title = element_text(size = 10))

pcoa_title3 = paste0(r,"Miseq Distance Matrix for  Sequence lane")
pcoa3 <- plot_ordination(glom, rt.pcoa,color="run_lane") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=run_lane)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title3)+theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 10))

beta_compartment_name <- paste0(r," Hiseq Compartment PCoA Beta Diversity.pdf")
ggsave(pcoa, filename=beta_compartment_name)

beta_regional_name <- paste0(r," Hiseq Expedition PCoA Beta Diversity.pdf")
ggsave(pcoa1, filename=beta_regional_name)

beta_family_name <- paste0(r," Hiseq Host Family PCoA Beta Diversity.pdf")
ggsave(pcoa2, filename=beta_family_name)
beta_run_name <- paste0(r," Hiseq Sequence lane PCoA Beta Diversity.pdf")
ggsave(pcoa3, filename=beta_run_name)

}

print(paste("Starting Loop across compartments"))
## Start for loop across compartments

compartment <- c("M","T","S")

for (i in compartment){
# Subset only corals from database
sub <-subset_samples(glom, tissue_compartment == print(paste(i)))
sub

paste(print("Creating glom mapping file"))
phyloseq::sample_data(sub)%>%
  as.data.frame() -> glom_mapping

#### create ASV tables by id ** This file will be used in microbe_neutral_compartment.R and picrust2_neutral_table_generator.R
print(paste("Generating Agglomerated ASV Table dataset..."))
phyloseq::otu_table(sub)%>%
  as.data.frame()%>%
  rownames_to_column("id") -> glom_otu_table

## Output .tsv from the otu table file
glom_otu_name <- paste0(i,"_glom_table.tsv")
write.table(glom_otu_table, file =glom_otu_name ,sep = "\t",row.names = FALSE)

#### create taxonomy tables by id  ** This 
paste(print("Printing Agglomerated Taxonomy Table"))

phyloseq::tax_table(sub)%>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  select(-c("Genus","Species"))-> glom_taxonomy

## Output .tsv from the taxonomy file
glom_taxonomy_name <- paste0(i,"_glom_taxonomy.tsv")
write.table(glom_taxonomy, file =glom_taxonomy_name,sep = "\t", row.names = FALSE)

## Output .csv from the taxonomy file
glom_mapping_name <- paste0(i,"_glom_metadata.tsv")
write.table(glom_mapping, file =glom_mapping_name,sep = "\t", row.names = TRUE)

########### Calculate Faiths Pd ############
print(paste("Calculating Faith Pd for host and microbiome from **glom** dataset"))

## pull out sample data 
phyloseq::sample_data(glom) %>%  
  group_by(sample_name_backup, expedition_number,BiologicalMatter) %>%
 as.data.frame() %>%
  select(sample_name_backup, expedition_number,BiologicalMatter) -> microbial_faith

## Microbial Phyloseq Pd analysis 
estimate_pd(glom) %>%
  as.data.frame() -> microbial_faith_pd

microbial_faith$faith_pd <- microbial_faith_pd$PD
microbial_faith$faith_SR <- microbial_faith_pd$SR


## Create data frame from sample data 
phyloseq::sample_data(glom) %>%  
  group_by(expedition_number, BiologicalMatter,Huang_Roy_tree_name) %>%
  as.data.frame() %>%
  select(expedition_number, BiologicalMatter,Huang_Roy_tree_name)-> test_df

## Create a new column titled eco to join expedition and biological matter
test_df$eco <-paste(test_df$expedition_number,test_df$BiologicalMatter, sep = "_")

## group and count total for each unique group maintaining NA values
test_df %>% group_by(eco, Huang_Roy_tree_name) %>%  summarise(counts=n()) %>%
  ungroup %>%
  complete(nesting(eco),
           nesting(Huang_Roy_tree_name),
          fill = list(quantity = 0)) -> test_table
# Fill NA with 0
test_table[is.na(test_table)] <-0

## Build Matrix
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

write.table(microbial_faith, file =paste0(i,"_","microbial_faithpd_table.txt"), sep = "\t",row.names =TRUE, col.names = TRUE)

write.table(coral_faithpd_reorded, file =paste0(i,"_","Host_faithpd_table.txt"), sep = "\t",row.names =TRUE, col.names = TRUE)
}

print(paste("Finished!"))

#