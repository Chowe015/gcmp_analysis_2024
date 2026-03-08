#Loading Required libraries
library(BiocManager)
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(microbiome)
library(stats4)
library(DESeq2)
library(vegan)
library(rstatix)


## Step 4 compare combined psudo tables to make comaprison between above, below and neutral taxa
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]
biosample <-args[5]


sink_name = paste0("./pseudo/",biosample,"/",biosample,"_Neutral_community_comparative_analysis_log.txt")
sink(sink_name,append=FALSE,split=TRUE)
####Import from .qza file into a phyloseq object

print(paste("feature_table",feature_table_path))
asv <- qza_to_phyloseq(features = feature_table_path)

#### Import Metadata read.table
metadata <- read.table(file = metadata_path,header=T, comment.char="",row.names = 1, sep="\t")

#### Import Tree file from biom output tree.nwk

print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)

#### Import taxonomy from biom output as .tsv format using read.table
print(paste("Loading Taxonomy text files from path:", taxonomy_path))
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

#biosample <-"mucus"
#asv <- qza_to_phyloseq(feature="./pseudo/mucus/mucus_feature_table.qza")
#metadata <-read.table(file="./pseudo/mucus/M_combined_metadata.txt",header=T, comment.char="",row.names = 1, sep="\t")
#tree <- read_tree("./tree.nwk")
#taxonomy <- read.table(file="./pseudo/mucus/mucus_pseudo_tax.tsv", sep = "\t", header = T ,row.names = 1 )
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

#print(paste("Loading metadata files from path:", metadata_path))
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
table(tax_table(phylo)[,"Kingdom"])
phylo

## Alpha Diversity Comparative Analysis and Boxplots
## Asssign metadata
sample_data(phylo) %>% as.data.frame() -> sample_meta

## Estimate Richness 
estimate_richness(phylo, measures= c("Observed","Shannon")) %>% as.data.frame() -> alpha
alpha$sample_name_backup <- sample_meta$sample_name_backup
alpha$neutrality <- sample_meta$neutrality
alpha$expedition_number <- sample_meta$expedition_number

#kruskal Wallis Test
krust_test <- kruskal.test(alpha$Observed ~ neutrality,data=alpha )
krust_test

#Wilcox pairwise test
pairwise_result <- pairwise.wilcox.test(alpha$Observed,alpha$neutrality, p.adjust.method="BH")
(pairwise_result)

## subset alpha results table
alpha %>% filter(neutrality=="above") %>% as.data.frame () -> above_alpha
alpha %>% filter(neutrality=="below") %>% as.data.frame () -> below_alpha
alpha %>% filter(neutrality=="neutral") %>% as.data.frame () -> neutral_alpha

## Print the median observed abundance
print(paste(biosample,"Median above microbes ==",median(above_alpha$Observed)))
print(paste(biosample,"Median below microbes ==",median(below_alpha$Observed)))
print(paste(biosample,"Median neutral microbes ==",median(neutral_alpha$Observed)))


## Plot alpha diversity Plot
print(paste("plotting richness alpha diversity scores"))
compare = list(c("below","above"),c("below","neutral"),c("above","neutral"))

alpha_plot_name <-paste0(biosample," Alpha Diversity Boxplot")

alpha_plot <- plot_richness(phylo, x="neutrality", color = "neutrality", title = alpha_plot_name, measures = ("Observed"),)+ geom_boxplot() +  theme_classic() + 
scale_color_manual(values=c("#E41A1C","#FF7F00","#999999")) + coord_cartesian(ylim=c(0,150))+
theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x.bottom = element_text(angle = -90)) +
stat_compare_means(method= "wilcox.test",comparisons = compare, symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))

## Print Neutral Model Plot
print(paste("Creating Neutral alpha diversity Box Plots"))
alpha_name <-paste0("./pseudo/",biosample,"/",biosample,"_Alpha_Diversity_Boxplot.pdf")
ggsave(alpha_plot, filename=alpha_name)

##### Differential Abundance Analysis ##########
print(paste("Running Differential Abundance Analysis..."))

## Assign comparison between treatments
compare = list(c("below","above"),c("below","neutral"),c("neutral","above"))

## Run for loop across each comparison
for (a in compare) {
  print(a)
  ## Subset samples for variables
  phylo.sub <- subset_samples(phylo,neutrality %in% a)
  
  
  ## Create a deseq-formatted matrix
  phylo.des <- phyloseq_to_deseq2(phylo.sub, ~neutrality)
  
  ## Caluclate geometric means prior to estimate size factor
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(phylo.des), 1, gm_mean)
  
  #Estimate size factor using geometric means calculation
  phylo.des = estimateSizeFactors(phylo.des, geoMeans = geoMeans)
  paste(a[1])
  ## Set the reference level * control 
  phylo.des$neutrality <- relevel(phylo.des$neutrality, ref = paste(a[1]))
  
  
  ##Run Deseq Command
  phylo.ds = DESeq(phylo.des)
  
  # Set the alpha thresholdd
  alpha = 0.05
  
  ## Evaluate Results
  res <- results(phylo.ds)
  
  ##order results by the smallest p-value
  resOrdered<-res[order(res$padj, res$baseMean),] 
  
  ##how many adjusted p-values < 0.05
  summary(resOrdered)
  l=sum(resOrdered$padj< 0.05, na.rm=TRUE) 
  
  ##get differentially expressed genes
  sig_res <- resOrdered[!is.na(resOrdered$padj),]
  sig_res = sig_res[sig_res$padj < 0.05,]
  
 ## Assign Taxonomy to DE results table
  sig_res %>% as.data.frame() -> sig_table
  sig_table$asv <- row.names(sig_table)
  tax_table(phylo) %>% as.data.frame()-> tax_list
  tax_list$asv <- row.names(tax_list)
  
  ## Merge significant table and taxonomy 
  inner_join(sig_table,tax_list, by = "asv") -> res_tab
  
  ## output Taxonomic results table
  deseq_tax_name <- paste0("./pseudo/",biosample,"/",biosample,"_",a[1],"_",a[2],"_taxonomic_deseq_results.tsv")
  (write.table(res_tab, file =deseq_tax_name,sep="\t", row.names = FALSE))
  
  # Create a table contains most significant features
  taxa_sig = rownames(sig_res[1:l, ])
  
  many_col <- c("#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
                "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8",
                "#4DAF4A" ,"#FF7F00","#984EA3", "#FFFF33", "#A65628","#1B9E77", "#F781BF" ,"#D53E4F","#F46D43", "#FDAE61", "#FEE08B",
                "#E6F598","#ABDDA4","#3288BD","cyan1","#66C2A5","#E69F00", "#56B4E9","#999999","#F0E442","#0072B2","#D55E00","red3",
                "yellow","skyblue","deeppink","cyan3","coral4","darkkhaki","brown1","#A6CEE3","chocolate","darkorchid2","#FF7F00",
                "#66A61E","#FDBF6F", "#FB9A99", "#CAB2D6",  "#B15928","#1B9E77","#D95F02", "#FFFF99","#E7298A", "#7570B3","#666666",
                "#E6AB02", "brown3","#A6761D", "#377EB8", "#4DAF4A" ,"#984EA3","#FFFF33","#FF7F00","black","#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
                "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8",
                "#4DAF4A" ,"#FF7F00","#984EA3", "#FFFF33", "#A65628","#1B9E77", "#F781BF" ,"#D53E4F","#F46D43", "#FDAE61", "#FEE08B","darkorchid2","#FF7F00",
                "#66A61E","#FDBF6F", "#FB9A99", "#CAB2D6",  "#B15928","#1B9E77","#D95F02", "#FFFF99","#E7298A", "#7570B3","#666666",
                "#E6AB02", "brown3","#A6761D", "#377EB8", "#4DAF4A" ,"#984EA3","#FFFF33","#FF7F00","black","#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
                "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8","#FF7F00","black","#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
                "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8")
  
  # Set to Order
  x = tapply(res_tab$log2FoldChange, res_tab$Class, function(x) max(x))
  x = sort(x,decreasing = TRUE)
  res_tab$Class = factor(as.character(res_tab$Class), levels=names(x))
  
  # Set to Family
  #x = tapply(res_tab$log2FoldChange, res_tab$Family, function(x) max(x))
  #x = sort(x,decreasing = TRUE)
  #res_tab$Family = factor(as.character(res_tab$Family), levels=names(x))
  
  # Deseq barplot 
  de_title <- paste0(biosample," ",a[1]," ",a[2]," ","Differential Abundance Plot")
  #DE_ofra_bar <-ggplot(data=res_tab, aes(x=Class, y=log2FoldChange, fill = Class))  +
  #geom_bar(stat="identity", position = "stack") + 
  #  geom_hline(yintercept = 0.0, color = "black", linewidth = 1)+ scale_fill_manual(values = c(many_col))+ labs( x=" Class", y="log2FoldChange") #+ coord_flip() + ggtitle(de_title)+ theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=6,face="bold"),axis.title=element_text(size=12))+ theme(legend.position="none")
  
  ## Create differential Expressions Dot Plot for OFRA
  de_neutral_plot <- ggplot(res_tab, aes(y=Class, x=log2FoldChange, color=Class)) + 
    geom_vline(xintercept = 0.0, color = "black", linewidth = 0.8)+ scale_color_manual(values = c(many_col)) + geom_point(size=3) + theme(legend.position="none") + ggtitle(de_title)+ theme(plot.title = element_text(hjust = 0.5))
  
  de_title_name <- paste0("./pseudo/",biosample,"/",biosample,"_",a[1],"_",a[2],"_DE_dotPlot.pdf")
  ggsave(de_neutral_plot,filename = de_title_name)
}


print("finished!")