#Loading Required libraries
library(BiocManager)
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(DESeq2)
library(ComplexHeatmap)
library(ade4)
library(vegan)
library(ggplot2)
library(microbiome)
library(dendextend)
library(RVAideMemoire)
library(ape)
library(RColorBrewer)
library(ggsignif)
library(ANCOMBC)
library(biomformat)
library(cowplot) 


#Step 1: Import qiime2 tables, mapping, taxonomy, tree and coral phylogeny. 
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
coral_tree_path <-args[5]

sink_name = paste0("gcmp_MiSeq_alpha_beta_diversity_log.txt")
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
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names=1, sep="\t")

#### Import Tree file from biom output .nwk

#print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)

coral_tree <-read_tree(coral_tree_path)

#### Import taxonomy from biom output as .tsv format using read.table
print(paste("Loading Taxonomy text files from path:", taxonomy_path))
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

##code referenced from Yan Hui: email me@yanh.org github: yanhui09**#####

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
  phyloseq::subset_taxa(!Family %in% c("Mitochondria","Chloroplast","mitochondria"))

# verify present taxa in table
table(tax_table(phylo_noking)[,"Kingdom"])

# Remove any samples with less than 1 read
phylo = prune_samples(sample_sums(phylo_noking)>1, phylo_noking)
phylo

### rarefy to even depth
print(paste("Generating rarefied Coral dataset..."))
rarefied = rarefy_even_depth(phylo, rngseed=111, sample.size=1000, replace=F, trimOTUs = TRUE)
print(rarefied)

#### Agglomerate taxa to family 
print(paste("Agglomerated Taxonomy to the Family Level"))
glom <- tax_glom(rarefied, taxrank = 'Family', NArm = TRUE)
print(glom)

# Inital Subset for agglomeration and remove outgroup or doubled samples
subject <- subset_samples(glom, outgroup=="n" & filter_col=="no")
subject

############## Analyze alpha diversity ####################################
print(paste("Calculating Alpha Diversity Results"))
## Generate Alpha Diversity data by estimating richness. *Note some richness estimates can't be calculated so call specific measures
gcmp_rich = estimate_richness(subject, split=T,measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))

# Subset metadata to assign to samples
phyloseq::sample_data(subject) %>% as.data.frame() -> sample_rich

# Assign metadata to alpha diversity results table
gcmp_rich$family <- sample_rich$family
gcmp_rich$BiologicalMatter <- sample_rich$BiologicalMatter
gcmp_rich$expedition_number <- sample_rich$expedition_number

# Set Factor for Compartments AND order by mucus tissue and skeleton
comp_level <- c("Coral Mucus","Coral Tissue","Coral Skeleton")
gcmp_rich$BiologicalMatter <- factor(gcmp_rich$BiologicalMatter, levels=comp_level)

# linear regression analysis with alpda diversity
## Multiple Analysis Across family+BiologicalMatter+expedition_number
alpha_lm_all <- lm(Shannon ~ family+BiologicalMatter+expedition_number, gcmp_rich)
summary(alpha_lm_all)

## Multiple regression Analysis Across Biology interactions
print(paste("Analyzing alpha diversity family and compartment regression interactions"))
alpha_lm_host <- lm(Shannon ~ family*BiologicalMatter, gcmp_rich)
summary(alpha_lm_host)

## ## Multiple regression Analysis Across Biology
print(paste("Analyzing alpha diversity family and compartment regression"))
alpha_lm_host2 <- lm(Shannon ~ family+BiologicalMatter, gcmp_rich)
summary(alpha_lm_host2)


## Linear regression with expedition_number
print(paste("Analyzing alpha diversity regression interactions between family and region"))
alpha_lm_exp_fam <- lm(Shannon ~ expedition_number*family, gcmp_rich)
summary(alpha_lm_exp_fam)

## Linear regression with expedition_number and family
print(paste("Analyzing alpha diversity regression between family and region"))
alpha_lm_exp_fam2 <- lm(Shannon ~ expedition_number+family, gcmp_rich)
summary(alpha_lm_exp_fam2)

## Linear regression with expedition_number
alpha_lm_exp <- lm(Shannon ~ expedition_number, gcmp_rich)
summary(alpha_lm_exp)

##Linear regression with Family
alpha_lm_fam <- lm(Shannon ~ family, gcmp_rich)
summary(alpha_lm_fam)

#Linear regression with Compartment mucus, tissue and skeleton
alpha_lm_comp <- lm(Shannon ~ BiologicalMatter, gcmp_rich)
summary(alpha_lm_comp)

## Wilcox pariwise comparison
comp.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, gcmp_rich$BiologicalMatter, p.adjust.method = "BH")
comp.shannon

region.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, gcmp_rich$expedition_number, p.adjust.method = "BH")
region.shannon

family.shannon <- pairwise.wilcox.test(gcmp_rich$Shannon, gcmp_rich$family, p.adjust.method = "BH")
family.shannon


## Alpha Diversity Box Plots: Biological Matter
alpha_compartment_all <- ggplot(gcmp_rich, aes(x=BiologicalMatter, y= Shannon, color=BiologicalMatter))+
  geom_boxplot(outlier.shape = NA,width=0.5)+ stat_boxplot(geom = "errorbar", width = 0.25) + scale_color_manual(values=c("#009197","#EC8C5C","#8E3B97"))+ geom_point(na.rm = TRUE) +
  theme_classic() + theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") +
  geom_signif(comparisons = list(c("Coral Mucus", "Coral Tissue"), c("Coral Skeleton","Coral Mucus"), c("Coral Tissue", "Coral Skeleton")), y_position = c(6,5,5.5))

## Alpha Diversity Box Plots: Coral Family
alpha_family_all <- ggplot(gcmp_rich, aes(x=family, y= Shannon, color=family)) +
  geom_boxplot(outlier.shape = NA,width=0.5)+ stat_boxplot(geom = "errorbar", width = 0.25) + scale_color_manual(values=many_col)+ geom_point(na.rm = TRUE) + theme_classic()+
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") 

## Alpha Diversity Box Plots: Sampling Expedition
alpha_region_all <- ggplot(gcmp_rich, aes(x=expedition_number, y= Shannon, color=expedition_number)) +
  geom_boxplot(outlier.shape = NA,width=0.5)+ stat_boxplot(geom = "errorbar", width = 0.25) + scale_color_manual(values=many_col)+ geom_point(na.rm = TRUE) + theme_classic()+
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + guides(col="none") 

## Output Alpha Diversity Plot plots as Pdf.
print(paste("Plotting Alpha Diversity Results"))
alpha_compartment_name <- paste0("Miseq Compartment Alpha Diversity.pdf")
ggsave(alpha_compartment_all, filename=alpha_compartment_name)

alpha_region_name <- paste0("Miseq GCMP Expedition Alpha Diversity.pdf")
ggsave(alpha_region_all, filename=alpha_region_name)

alpha_family_name <- paste0("Miseq Host Family Alpha Diversity.pdf")
ggsave(alpha_family_all, filename=alpha_family_name)


print(paste("Calculating Bray-Curtis Matrices "))
# All samples
bray_curtis = phyloseq::distance(subject, method="bray")
print(paste("Calculating Unifrac Matrices "))
unifrac = phyloseq::distance(subject, method="unifrac")
print(paste("Calculating Weighted Matrices "))
weighted_unifrac = phyloseq::distance(subject, method="wUnifrac")

# Create data frame for downstream analysis
all_df <- data.frame(sample_data(subject))

############# PERMANOVA analysis###################
print(paste("Calculating Beta Diversity Results"))

distance_methods <-c("bray_curtis","unifrac","weighted_unifrac")

set.seed(129)
# This is what is called a for loop.        
for (i in distance_methods){ 
  form <- as.formula(paste(i,"BiologicalMatter", sep="~"))
  print(form)
 adonis2(form, data=all_df)-> result
  print(result)
}

set.seed(129)
# This is what is called a for loop.        
for (i in distance_methods){ 
  form <- as.formula(paste(i,"expedition_number", sep="~"))
  print(form)
 adonis2(form, data=all_df)-> result
  print(result)
}

set.seed(129)
# This is what is called a for loop.        
for (i in distance_methods){ 
  form <- as.formula(paste(i,"family", sep="~"))
  print(form)
  adonis2(form, data=all_df)-> result
  print(result)
}

print(paste("Plotting Beta Diversity Results"))

###Calculate the PCoA on Compartment, Region and Family 
for (r in distance_methods){
print(paste0(r,"Ordinating Distance Matrix"))

# Ordinate distances
rt.pcoa = ordinate(subject, method="PCoA", distance=r)

pcoa_title = paste0(r,"Miseq Distance Matrix for  Compartments")
## Plot PcoA for Compartment
pcoa <- plot_ordination(subject, rt.pcoa,color="BiologicalMatter") + scale_color_manual(values=c("#009197","#8E3B97","#EC8C5C"))+ geom_point(size=3) + theme_bw() +
stat_ellipse(level = 0.95,type = "norm",aes(group=BiologicalMatter)) + theme(legend.text=element_text(size=12)) +
theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title)+theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size = 10))
## Plot PcoA for Expedition
pcoa_title1 = paste0(r,"Miseq Distance Matrix for  Expeditions")
pcoa1 <- plot_ordination(subject, rt.pcoa,color="expedition_number") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=expedition_number)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title1)+theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size = 10))
## Plot PcoA for Family
pcoa_title2 = paste0(r,"Miseq Distance Matrix for  Host Family")
pcoa2 <- plot_ordination(subject, rt.pcoa,color="family") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
 stat_ellipse(level = 0.95,type = "norm",aes(group=family)) + theme(legend.text=element_text(size=12)) +
 theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title2)+theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size = 10))
## Plot PcoA for Sequence Lane
pcoa_title3 = paste0(r,"Miseq Distance Matrix for  Sequence lane")
pcoa3 <- plot_ordination(subject, rt.pcoa,color="run_lane") + scale_color_manual(values=many_col)+ geom_point(size=3) + theme_bw() +
  stat_ellipse(level = 0.95,type = "norm",aes(group=run_lane)) + theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14,face="bold") +
          theme(text = element_text(size = 20,face = "bold")))  + theme(text = element_text(size = 20,face = "bold"))+
  ggtitle(pcoa_title3)+theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 10))

## Output Beta Diversity Plots as Pdf.
beta_compartment_name <- paste0(r,"  Miseq Compartment PCoA Beta Diversity.pdf")
ggsave(pcoa, filename=beta_compartment_name)

beta_regional_name <- paste0(r," Miseq Expedition PCoA Beta Diversity.pdf")
ggsave(pcoa1, filename=beta_regional_name)

beta_family_name <- paste0(r," Miseq Host Family PCoA Beta Diversity.pdf")
ggsave(pcoa2, filename=beta_family_name)

beta_run_name <- paste0(r," Miseq sequence lane PCoA Beta Diversity.pdf")
ggsave(pcoa3, filename=beta_run_name)


}

print(paste("Finished!"))
