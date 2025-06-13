#read in contributing otus table
otus_contributing <- read.table("/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Potentially_Pathogenic/Mucus_neutral_model_below_table-97_Potentially_Pathogenic/otu_contributions/contributing_otus.txt")

#read in gg_taxonomy table
gg_taxonomy <- read.table("/Users/yifanli/Documents/PhD/Projects/GCMP/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt", sep="\t", row.names=1, check=F, quote='')

#read in normalized otu table for Aerobic above
otu_table <- read.table("/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Potentially_Pathogenic/Mucus_neutral_model_below_table-97_Potentially_Pathogenic/normalized_otus/16s_normalized_otus.txt", row.names=1, check=F, quote='')

#keep only otus in the gg taxonomy table that are in the otu table
#these tables will be in the same otu order
#OTUs are rows
otus_keep <- intersect(rownames(otu_table),rownames(gg_taxonomy))
gg_taxonomy <- gg_taxonomy[otus_keep,, drop=F]

#multiple the otu_table by the boolean table for otus contributing to 
# the trait

positive_otus <- intersect(rownames(otu_table),
                           rownames(otus_contributing)[which(otus_contributing$Potentially_Pathogenic=="TRUE")])
positive_otu_table <- otu_table[rownames(otu_table) %in% positive_otus,]

#subset the gg taxonomy to the level specified (family=5)
taxa_level=5
names_split <- array(dim=c(length(gg_taxonomy[,1]), 7))
rownames(names_split) <- otus_keep
otu_names <- as.character(gg_taxonomy[,1])
for(i in 1:length(otu_names)){
  names_split[i,] <- strsplit(otu_names[i], ";", fixed=T)[[1]]
}
otu_names <- names_split[,taxa_level]
for(i in 1:length(otu_names)){
  otu_names[i] <- strsplit(otu_names[i], "__", fixed=T)[[1]][2]
}
names_split[,taxa_level] <- otu_names
for(i in 1:nrow(names_split)){
  if(is.na(names_split[i,taxa_level])){
    if(taxa_level > 1){
      names_split[i, taxa_level] <- names_split[i, taxa_level -1]
    } else {
      names_split[i, taxa_level] <- "unknown"
    }
  }
}

#aggregate to the same taxa
#you must t() again, to have samples as columns
positive_otu_table <- t(sapply(by(positive_otu_table,rownames(positive_otu_table),colSums),identity))

#add taxonomy as the rownames in the otu table
positive_taxa <- intersect(rownames(names_split),
                           rownames(positive_otu_table))
names_split <- names_split[rownames(names_split) %in% positive_taxa,]
rownames(positive_otu_table) <- names_split[taxa_level]

#set colors for plotting and legend creation
#get palette 1 from R ColorBrewer
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9,'Set1'))

#create as many colors as there are taxa, and name them with the taxa names
cols2 <- cols(length(rownames(otu_table)))
names(cols2) <- unique(rownames(otu_table))
cols2 <- c(cols2,"#C0C0C0")
names(cols2)[length(cols2)] <- "Other"

taxa_list <- c()

xlabel <- "Relative Abundance"


#melt the otu_table and collapse by Taxa
library(reshape2)
library(plyr)
melted_otu_table <- melt(positive_otu_table)
colnames(melted_otu_table) <- c("Taxa", "SampleID", "Count")
melted_otu_table$SampleID <- as.factor(1)
taxa_collapsed <- ddply(melted_otu_table, .(Taxa, SampleID),summarize, 
                        Count = mean(Count))
taxa_collapsed$Taxa <- as.character(taxa_collapsed$Taxa)
taxa_collapsed$Count <- taxa_collapsed$Count/length(otu_table)
taxa_collapsed[2,] <- c('Non-contributing_OTUs', 1, 1-taxa_collapsed$Count)
class(taxa_collapsed$Count) <- "numeric"

#make the plot
library(ggplot2)
taxa_plot <- NULL
taxa_plot <- ggplot(taxa_collapsed, aes(fill=Taxa, y=Count, x=SampleID)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()
taxa_plot
