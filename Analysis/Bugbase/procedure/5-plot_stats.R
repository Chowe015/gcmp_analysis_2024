# Import csv tables of predicted phenotypes of a particular trait.
library(reshape2)

## Print commands to command line
options(echo=TRUE)

## Get user input and assign to variables
cat("Please enter the trait of interest.")
trait <- readLines(con = "stdin", n = 1)

wd = paste("/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/", trait, "/", sep="")
setwd(wd)

### Mucus
mucus_above_tsv <- paste(wd, "Mucus_Neutralmodel_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

mucus_neutral_tsv <- paste(wd, "Mucus_Neutralmodel_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

mucus_below_tsv <- paste(wd, "Mucus_Neutralmodel_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

mucus_all_tsv <- paste(wd, "Mucus_Neutralmodel_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

#mucus_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#mucus_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#mucus_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#mucus_all_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_table-97_Aerobic/predicted_phenotypes/predictions.txt"

M_above = read.delim(file=mucus_above_tsv, header = T, sep="\t")
colnames(M_above) <- c("sampleID","M_above")
M_above = melt(M_above)
M_above$compartment = "Mucus"
M_neutral = read.delim(file=mucus_neutral_tsv, header = T, sep="\t")
colnames(M_neutral) <- c("sampleID","M_neutral")
M_neutral = melt(M_neutral)
M_neutral$compartment = "Mucus"
M_below = read.delim(file=mucus_below_tsv, header = T, sep="\t")
colnames(M_below) <- c("sampleID","M_below")
M_below = melt(M_below)
M_below$compartment = "Mucus"

M_all = read.delim(file=mucus_all_tsv, header = T, sep="\t")
colnames(M_all) <- c("sampleID","M_all")
M_all = melt(M_all)
M_all$compartment = "Mucus"

### Tissue
tissue_above_tsv <- paste(wd, "Tissue_Neutralmodel_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

tissue_neutral_tsv <- paste(wd, "Tissue_Neutralmodel_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

tissue_below_tsv <- paste(wd, "Tissue_Neutralmodel_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

tissue_all_tsv <- paste(wd, "Tissue_Neutralmodel_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

#tissue_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#tissue_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#tissue_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#tissue_all_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_table-97_Aerobic/predicted_phenotypes/predictions.txt"

T_above = read.delim(file=tissue_above_tsv, header = T, sep="\t")
colnames(T_above) <- c("sampleID","T_above")
T_above = melt(T_above)
T_above$compartment = "Tissue"
T_neutral = read.delim(file=tissue_neutral_tsv, header = T, sep="\t")
colnames(T_neutral) <- c("sampleID","T_neutral")
T_neutral = melt(T_neutral)
T_neutral$compartment = "Tissue"
T_below = read.delim(file=tissue_below_tsv, header = T, sep="\t")
colnames(T_below) <- c("sampleID","T_below")
T_below = melt(T_below)
T_below$compartment = "Tissue"

T_all = read.delim(file=tissue_all_tsv, header = T, sep="\t")
colnames(T_all) <- c("sampleID","T_all")
T_all = melt(T_all)
T_all$compartment = "Tissue"

### Skeleton
skeleton_above_tsv <- paste(wd, "Skeleton_Neutralmodel_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

skeleton_neutral_tsv <- paste(wd, "Skeleton_Neutralmodel_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

skeleton_below_tsv <- paste(wd, "Skeleton_Neutralmodel_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

skeleton_all_tsv <- paste(wd, "Skeleton_Neutralmodel_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")

#skeleton_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#skeleton_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#skeleton_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#skeleton_all_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_table-97_Aerobic/predicted_phenotypes/predictions.txt"

S_above = read.csv(file=skeleton_above_tsv, header = T, sep="\t")
colnames(S_above) <- c("sampleID","S_above")
S_above = melt(S_above)
S_above$compartment = "Skeleton"
S_neutral = read.csv(file=skeleton_neutral_tsv, header = T, sep="\t")
colnames(S_neutral) <- c("sampleID","S_neutral")
S_neutral = melt(S_neutral)
S_neutral$compartment = "Skeleton"
S_below = read.csv(file=skeleton_below_tsv, header = T, sep="\t")
colnames(S_below) <- c("sampleID","S_below")
S_below = melt(S_below)
S_below$compartment = "Skeleton"

S_all = read.delim(file=skeleton_all_tsv, header = T, sep="\t")
colnames(S_all) <- c("sampleID","S_all")
S_all = melt(S_all)
S_all$compartment = "Skeleton"


### Plot figures for pooled MST ###

# Combine all predicted phenotypes into one table.
library(dplyr)

#merge all data frames in list
all_aerobic_pooled <- rbind(M_all, T_all, S_all)
all_aerobic_pooled$sampleID <- substr(all_aerobic_pooled$sampleID, 1, nchar(all_aerobic_pooled$sampleID)-2)

# Retain only samples that are shared between all compartments and model fits.
M_all$sampleID <- substr(M_all$sampleID, 1, nchar(M_all$sampleID)-2)
T_all$sampleID <- substr(T_all$sampleID, 1, nchar(T_all$sampleID)-2)
S_all$sampleID <- substr(S_all$sampleID, 1, nchar(S_all$sampleID)-2)
shared_all_list=Reduce(intersect, list(M_all$sampleID,
                                  T_all$sampleID,
                                  S_all$sampleID))
shared_all <- filter(all_aerobic_pooled,
                       sampleID %in% shared_all_list)

# Plot violin.
library(ggplot2)

## Basic violin plot
p <- ggplot(shared_all, aes(x=variable, y=value)) + theme_classic() +
  geom_violin(aes(fill=compartment), trim=F) +
  scale_fill_manual(breaks = c("Mucus", "Tissue", "Skeleton"), values=c("#009197","#EC8C5C","#8E3B97"))
p

## Add stats (Kruskal-Wallis)
library(rstatix)
library(ggpubr)

# kw_results = all_aerobic_pooled %>% kruskal_test(value ~ variable)

# Wilcoxon pair-wise comparisons
### Remove ties from Wilcoxon stats ###
# Step 1: Store dataframes in a named list
df_list <- list(
  M_all = M_all,
  T_all = T_all,
  S_all = S_all)


# List of category names
categories <- names(df_list)

# Step 2: Double for loop to compare each pair of categories
cols = c("comparison", "p_value")
wilcoxon_pooled = data.frame(matrix(ncol = 2, nrow = 0))
colnames(wilcoxon_pooled) = cols

for (i in 1:length(categories)) {
  for (j in 1:length(categories)) {
    
    
    #Skip self comparisons
    if (i == j) next
    
    
    cat1_name <- categories[i]
    cat2_name <- categories[j]
    
    
    # Extract the dataframes from the list of dataframes
    # based on the current categories
    df1 <- df_list[[cat1_name]]
    df2 <- df_list[[cat2_name]]
    
    # Print which categories are being compared
    cat("Comparing", categories[i], "vs.", categories[j], "\n")
    
    
    
    
    # Optional: Check the first few rows of each dataframe (for demonstration)
    print(head(df1))
    print(head(df2))
    
    
    
    merged_data <- merge(df1, df2, by = "sampleID")
    print(merged_data)
    
    # --- Your comparison code here (e.g., Wilcoxon test, merging, etc.) ---
    
    if (nrow(merged_data) > 0) {
      non_equal_pairs <- merged_data$value.x != merged_data$value.y
      print(non_equal_pairs)
      
      wilcox_result <- wilcox.test(merged_data$value.x[non_equal_pairs],
                                   merged_data$value.y[non_equal_pairs],
                                   paired=TRUE)
      print(cat("Comparing", categories[i], "vs.", categories[j], "\n"))
      print(wilcox_result)
      
      dir.create(paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait,sep=""))
      outputfile=paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait, "/", trait,"_wilcoxon.txt", sep="")
      
      capture.output(print(cat("Comparing", categories[i], "vs.", categories[j], "\n")),
                     file = outputfile,
                     append=T)
      capture.output(wilcox_result, 
                     file = outputfile,
                     append=T)
      
      comparison = paste(categories[i], "vs.", categories[j])
      p_value = print(wilcox_result$p.value)
      wilcoxon_pooled = rbind(wilcoxon_pooled, c(comparison, p_value))
      colnames(wilcoxon_pooled) = cols
    }
  }
}

# Remove duplicate comparisons
keep=c("M_all vs. T_all", "M_all vs. S_all", "S_all vs. T_all")
wilcoxon_pooled = subset(wilcoxon_pooled, wilcoxon_pooled$comparison == keep)

# Perform FDR on p-values
wilcoxon_pooled$p.adjusted = signif(p.adjust(wilcoxon_pooled$p_value, method = "fdr"), 5)

# Split comparisons into two columns
wilcoxon_pooled=separate_wider_delim(wilcoxon_pooled, comparison, " ", names=c("compare1", "vs", "compare2"))

p <- ggplot(shared_all, aes(x=variable, y=value)) + theme_classic() +
  geom_violin(aes(fill=compartment), trim=F) +
  scale_fill_manual(breaks = c("Mucus", "Tissue", "Skeleton"), values=c("#009197","#EC8C5C","#8E3B97"))
p = p + 
  geom_bracket(
    xmin = wilcoxon_pooled$compare1,
    xmax = wilcoxon_pooled$compare2,
    label = wilcoxon_pooled$p.adjusted,
    y.position = c(1.5, 1.8, 2.1)
  ) + 
  geom_jitter(shape=16, size=1, position=position_jitter(0.2)) +
  stat_summary(fun=median, geom="point", shape=23, size=3, fill="white", color="grey50") + stat_summary(fun=mean, geom="point", size=2, color ="white") +
  stat_compare_means(label.x = 2.5, label.y = 2.5)
p = p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", shape=21, color="grey50", fill="grey90", size =0.8) +
  labs(x = "Compartments, pooled all model fit", y = "Relative abundance of predicted phenotypes")
p

savefile=paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait, "/", trait,"_pooled.pdf", sep="")
ggsave(savefile, p, create.dir = TRUE,
       width = 12.6,
       height = 9.06,
       units = "in")

### Plot figures for split MST ###

# Combine all predicted phenotypes into one table.
library(dplyr)

# Merge all data frames in list
all_aerobic_split = rbind(M_above,M_neutral,M_below,T_above,T_neutral,T_below,S_above,S_neutral,S_below)
all_aerobic_split$sampleID <- substr(all_aerobic_split$sampleID, 1, nchar(all_aerobic_split$sampleID)-2)

# Plot violin.
library(ggplot2)

## Basic violin plot
p <- ggplot(all_aerobic_split, aes(x=variable, y=value)) + theme_classic() +
  geom_violin(aes(fill=compartment), trim=F) +
  scale_fill_manual(breaks = c("Mucus", "Tissue", "Skeleton"), values=c("#009197","#EC8C5C","#8E3B97"))
p

#savefile=paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait, "/", trait,"_split_nostats.pdf", sep="")
#ggsave(savefile, p, create.dir = TRUE,
#       width = 12.6,
#       height = 9.06,
#       units = "in")

## Add stats (Wilcoxon)
library(rstatix)
library(ggpubr)

#wilcoxon_results = all_aerobic_split %>% pairwise_wilcox_test(value ~ variable, p.adjust.method="fdr")
#wilcoxon_sig = wilcoxon_results %>% filter(p.adj.signif != "ns")

# Wilcoxon pair-wise comparisons
### Remove ties from Wilcoxon stats ###
# Step 1: Store dataframes in a named list
df_list <- list(
  M_above = M_above,
  M_neutral = M_neutral,
  M_below = M_below,
  T_above = T_above,
  T_neutral = T_neutral,
  T_below = T_below,
  S_above = S_above,
  S_neutral = S_neutral,
  S_below = S_below)


# List of category names
categories <- names(df_list)

# Step 2: Double for loop to compare each pair of categories
cols = c("comparison", "p_value")
wilcoxon_pooled = data.frame(matrix(ncol = 2, nrow = 0))
colnames(wilcoxon_pooled) = cols

for (i in 1:length(categories)) {
  for (j in 1:length(categories)) {
    
    
    #Skip self comparisons
    if (i == j) next
    
    
    cat1_name <- categories[i]
    cat2_name <- categories[j]
    
    
    # Extract the dataframes from the list of dataframes
    # based on the current categories
    df1 <- df_list[[cat1_name]]
    df2 <- df_list[[cat2_name]]
    
    # Print which categories are being compared
    cat("Comparing", categories[i], "vs.", categories[j], "\n")
    
    
    
    
    # Optional: Check the first few rows of each dataframe (for demonstration)
    print(head(df1))
    print(head(df2))
    
    
    
    merged_data <- merge(df1, df2, by = "sampleID")
    print(merged_data)
    
    # --- Your comparison code here (e.g., Wilcoxon test, merging, etc.) ---
    
    if (nrow(merged_data) > 0) {
      non_equal_pairs <- merged_data$value.x != merged_data$value.y
      print(non_equal_pairs)
      
      wilcox_result <- wilcox.test(merged_data$value.x[non_equal_pairs],
                                   merged_data$value.y[non_equal_pairs],
                                   paired=TRUE)
      print(cat("Comparing", categories[i], "vs.", categories[j], "\n"))
      print(wilcox_result)
      
      dir.create(paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait,sep=""))
      outputfile=paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait, "/", trait,"_wilcoxon.txt", sep="")
      
      capture.output(print(cat("Comparing", categories[i], "vs.", categories[j], "\n")),
                     file = outputfile,
                     append=T)
      capture.output(wilcox_result, 
                     file = outputfile,
                     append=T)
      
      comparison = paste(categories[i], "vs.", categories[j])
      p_value = print(wilcox_result$p.value)
      wilcoxon_pooled = rbind(wilcoxon_pooled, c(comparison, p_value))
      colnames(wilcoxon_pooled) = cols
    }
  }
}

# Remove duplicate comparisons
keep=c("M_above vs. M_neutral", 
       "M_above vs. M_below",  
       "M_above vs. T_above", 
       "M_above vs. T_neutral",  
       "M_above vs. T_below",  
       "M_above vs. S_above",  
       "M_above vs. S_neutral",  
       "M_above vs. S_below",  
       "M_neutral vs. M_below",  
       "M_neutral vs. T_above",  
       "M_neutral vs. T_neutral",  
       "M_neutral vs. T_below",  
       "M_neutral vs. S_above",  
       "M_neutral vs. S_neutral",  
       "M_neutral vs. S_below",  
       "M_below vs. T_above",  
       "M_below vs. T_neutral",  
       "M_below vs. T_below",  
       "M_below vs. S_above",  
       "M_below vs. S_neutral",  
       "M_below vs. S_below",  
       "T_above vs. T_neutral",  
       "T_above vs. T_below",  
       "T_above vs. S_above",  
       "T_above vs. S_neutral",  
       "T_above vs. S_below",  
       "T_neutral vs. T_below",  
       "T_neutral vs. S_above",  
       "T_neutral vs. S_neutral",  
       "T_neutral vs. S_below",  
       "T_below vs. S_above",  
       "T_below vs. S_neutral",  
       "T_below vs. S_below",  
       "S_above vs. S_neutral",  
       "S_above vs. S_below",  
       "S_neutral vs. S_below")
wilcoxon_pooled = subset(wilcoxon_pooled, wilcoxon_pooled$comparison == keep)

# Perform FDR on p-values
wilcoxon_pooled$p.adjusted = signif(p.adjust(wilcoxon_pooled$p_value, method = "fdr"), 5)

# Split comparisons into two columns
wilcoxon_pooled=separate_wider_delim(wilcoxon_pooled, comparison, " ", names=c("compare1", "vs", "compare2"))

p = p + stat_pvalue_manual(wilcoxon_sig, y.position=1.6, step.increase=0.07) + geom_jitter(shape=16, size=1, position=position_jitter(0.2)) +
  stat_summary(fun=median, geom="point", shape=23, size=3, fill="white", color="grey50") + stat_summary(fun=mean, geom="point", size=2, color ="white") 
p = p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", shape=21, color="grey50", fill="grey90", size =0.8) +
  labs(x = "Compartments and respective model fit", y = "Relative abundance of predicted phenotypes") +
  stat_compare_means(label.x = 7.8, label.y = 4.8)
p

savefile=paste("/Users/yifanli/Documents/PhD/Projects/GCMP/R_figures/", trait, "/", trait,"_split.pdf", sep="")
ggsave(savefile, p, create.dir = TRUE,
       width = 12.6,
       height = 9.06,
       units = "in")
