# Import csv tables of predicted phenotypes of a particular trait.
library(reshape2)

## Print commands to command line
options(echo=TRUE)

## Get user input and assign to variables
cat("Please enter the trait of interest.")
trait <- readLines(con = "stdin", n = 1)

### Mucus
cat("Please enter path to predicted phenotypes file of bacteria in mucus compartment above the neutral model line.")
mucus_above_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in mucus compartment on the neutral model line.")
mucus_neutral_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in mucus compartment below the neutral model line.")
mucus_below_tsv <- readLines(con = "stdin", n = 1)

#mucus_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#mucus_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#mucus_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Mucus_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"

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

### Tissue
cat("Please enter path to predicted phenotypes file of bacteria in tissue compartment above the neutral model line.")
tissue_above_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in tissue compartment on the neutral model line.")
tissue_neutral_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in tissue compartment below the neutral model line.")
tissue_below_tsv <- readLines(con = "stdin", n = 1)

#tissue_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#tissue_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#tissue_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Tissue_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"

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

### Skeleton
cat("Please enter path to predicted phenotypes file of bacteria in skeleton compartment above the neutral model line.")
skeleton_above_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in skeleton compartment on the neutral model line.")
skeleton_neutral_tsv <- readLines(con = "stdin", n = 1)

cat("Please enter path to predicted phenotypes file of bacteria in skeleton compartment below the neutral model line.")
skeleton_below_tsv <- readLines(con = "stdin", n = 1)

#skeleton_above_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_above_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#skeleton_neutral_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_neutral_table-97_Aerobic/predicted_phenotypes/predictions.txt"
#skeleton_below_tsv = "/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/Aerobic/Skeleton_Neutralmodel_below_table-97_Aerobic/predicted_phenotypes/predictions.txt"

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

# Combine all predicted phenotypes into one table.
library(dplyr)

#merge all data frames in list
all_aerobic <- rbind(M_above,M_neutral,M_below,T_above,T_neutral,T_below,S_above,S_neutral,S_below)

# Plot violin.
library(ggplot2)

## Basic violin plot
p <- ggplot(all_aerobic, aes(x=variable, y=value)) + theme_classic() +
  geom_violin(aes(fill=compartment), trim=F) +
  scale_fill_manual(values=c("#009197","#EC8C5C","#8E3B97"))
p

## Add stats (Wilcoxon)
library(rstatix)
library(ggpubr)

wilcoxon_results = all_aerobic %>% pairwise_wilcox_test(value ~ variable)
wilcoxon_sig = wilcoxon_results %>% filter(p.adj.signif != "ns")

p = p + stat_pvalue_manual(wilcoxon_sig, y.position=1.6, step.increase=0.1) + geom_jitter(shape=16, size=1, position=position_jitter(0.2)) +
  stat_summary(fun=median, geom="point", shape=23, size=3, fill="white", color="white") + stat_summary(fun=mean, geom="point", size=2, color ="white") 
p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="grey") +
  labs(x = "Compartments and respective model fit", y = "Relative abundance of predicted phenotypes")
