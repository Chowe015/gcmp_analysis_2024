# Z-score normalization and combine all datasets. #

library(reshape2)
library(tidyr)
library(dplyr)

all_traits = c("Aerobic", "Anaerobic","Contains_Mobile_Elements", "Facultatively_Anaerobic", "Forms_Biofilms", "Gram_Negative", "Gram_Positive", "Potentially_Pathogenic", "Stress_Tolerant", "M00175_Nitrogen_fixation_nitrogen_ammonia_", "M00176_Sulfur_reduction_sulfate_H2S_", "M00453_QseC_QseB_quorum_sensing_two_component_regulatory_system_", "M00506_CheA_CheYBV_chemotaxis_two_component_regulatory_system_", "M00513_LuxQN_CqsS_LuxU_LuxO_quorum_sensing_two_component_regulatory_system_")

for (i in all_traits){
  
    trait = i
    
    wd = paste("/Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/", trait, "/", sep="")
    setwd(wd)
    
    ### Mucus
    mucus_above_tsv <- paste(wd, "Mucus_neutral_model_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    mucus_neutral_tsv <- paste(wd, "Mucus_neutral_model_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    mucus_below_tsv <- paste(wd, "Mucus_neutral_model_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    mucus_all_tsv <- paste(wd, "Mucus_neutral_model_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    
    M_above = read.delim(file=mucus_above_tsv, header = T, sep="\t")
    colnames(M_above) <- c("sampleID","M_above")
    M_above = melt(M_above)
    M_above$compartment = "Mucus"
    M_above$model_fit = "Above"
    M_above$sampleID <- substr(M_above$sampleID, 1, nchar(M_above$sampleID)-2)
    M_above$trait = i
    
    M_neutral = read.delim(file=mucus_neutral_tsv, header = T, sep="\t")
    colnames(M_neutral) <- c("sampleID","M_neutral")
    M_neutral = melt(M_neutral)
    M_neutral$compartment = "Mucus"
    M_neutral$model_fit = "Neutral"
    M_neutral$sampleID <- substr(M_neutral$sampleID, 1, nchar(M_neutral$sampleID)-2)
    M_neutral$trait = i
    
    M_below = read.delim(file=mucus_below_tsv, header = T, sep="\t")
    colnames(M_below) <- c("sampleID","M_below")
    M_below = melt(M_below)
    M_below$compartment = "Mucus"
    M_below$model_fit = "Below"
    M_below$sampleID <- substr(M_below$sampleID, 1, nchar(M_below$sampleID)-2)
    M_below$trait = i
    
    ### Tissue
    tissue_above_tsv <- paste(wd, "Tissue_neutral_model_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    tissue_neutral_tsv <- paste(wd, "Tissue_neutral_model_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    tissue_below_tsv <- paste(wd, "Tissue_neutral_model_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    tissue_all_tsv <- paste(wd, "Tissue_neutral_model_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    
    T_above = read.delim(file=tissue_above_tsv, header = T, sep="\t")
    colnames(T_above) <- c("sampleID","T_above")
    T_above = melt(T_above)
    T_above$compartment = "Tissue"
    T_above$model_fit = "Above"
    T_above$sampleID <- substr(T_above$sampleID, 1, nchar(T_above$sampleID)-2)
    T_above$trait = i
    
    T_neutral = read.delim(file=tissue_neutral_tsv, header = T, sep="\t")
    colnames(T_neutral) <- c("sampleID","T_neutral")
    T_neutral = melt(T_neutral)
    T_neutral$compartment = "Tissue"
    T_neutral$model_fit = "Neutral"
    T_neutral$sampleID <- substr(T_neutral$sampleID, 1, nchar(T_neutral$sampleID)-2)
    T_neutral$trait = i
    
    T_below = read.delim(file=tissue_below_tsv, header = T, sep="\t")
    colnames(T_below) <- c("sampleID","T_below")
    T_below = melt(T_below)
    T_below$compartment = "Tissue"
    T_below$model_fit = "Below"
    T_below$sampleID <- substr(T_below$sampleID, 1, nchar(T_below$sampleID)-2)
    T_below$trait = i
    
    ### Skeleton
    skeleton_above_tsv <- paste(wd, "Skeleton_neutral_model_above_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    skeleton_neutral_tsv <- paste(wd, "Skeleton_neutral_model_neutral_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    skeleton_below_tsv <- paste(wd, "Skeleton_neutral_model_below_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    skeleton_all_tsv <- paste(wd, "Skeleton_neutral_model_table-97_", trait, "/predicted_phenotypes/predictions.txt", sep="")
    
    S_above = read.delim(file=skeleton_above_tsv, header = T, sep="\t")
    colnames(S_above) <- c("sampleID","S_above")
    S_above = melt(S_above)
    S_above$compartment = "Skeleton"
    S_above$model_fit = "Above"
    S_above$sampleID <- substr(S_above$sampleID, 1, nchar(S_above$sampleID)-2)
    S_above$trait = i
    
    S_neutral = read.delim(file=skeleton_neutral_tsv, header = T, sep="\t")
    colnames(S_neutral) <- c("sampleID","S_neutral")
    S_neutral = melt(S_neutral)
    S_neutral$compartment = "Skeleton"
    S_neutral$model_fit = "Neutral"
    S_neutral$sampleID <- substr(S_neutral$sampleID, 1, nchar(S_neutral$sampleID)-2)
    S_neutral$trait = i
    
    S_below = read.delim(file=skeleton_below_tsv, header = T, sep="\t")
    colnames(S_below) <- c("sampleID","S_below")
    S_below = melt(S_below)
    S_below$compartment = "Skeleton"
    S_below$model_fit = "Below"
    S_below$sampleID <- substr(S_below$sampleID, 1, nchar(S_below$sampleID)-2)
    S_below$trait = i
    
    all_comp_fit = rbind(M_above, M_neutral, M_below, T_above, T_neutral, T_below, S_above, S_neutral, S_below)
    all_comp_fit$Z = scale(all_comp_fit$value)
    
    assign(as.vector(paste0(i)), all_comp_fit)
}

all_traits_tables = rbind(Aerobic, Anaerobic, Contains_Mobile_Elements, Facultatively_Anaerobic, Forms_Biofilms, Gram_Negative, Gram_Positive, Potentially_Pathogenic, Stress_Tolerant, M00175_Nitrogen_fixation_nitrogen_ammonia_, M00176_Sulfur_reduction_sulfate_H2S_, M00453_QseC_QseB_quorum_sensing_two_component_regulatory_system_, M00506_CheA_CheYBV_chemotaxis_two_component_regulatory_system_, M00513_LuxQN_CqsS_LuxU_LuxO_quorum_sensing_two_component_regulatory_system_)

# Calculate TukeyHSD
library(multcompView)

mod1 = all_traits_tables %>% group_by(trait) %>% do(model = lm(Z ~ variable, data = .))
all_traits_tables$conditions = paste(all_traits_tables$trait, all_traits_tables$variable, sep="_")

for (i in 1:14){
              Letters1 <- multcompLetters4(mod1$model[[i]], TukeyHSD(aov(mod1$model[[i]])))
              cld = as.data.frame.list(Letters1$variable)
              
              all_traits_tables = all_traits_tables %>%
                mutate(TukeyHSD = case_when(
                  conditions == paste(mod1$trait[i], "S_below", sep="_") ~ cld[row.names(cld)=="S_below",]$Letters,
                  conditions == paste(mod1$trait[i], "S_neutral", sep="_") ~ cld[row.names(cld)=="S_neutral",]$Letters,
                  conditions == paste(mod1$trait[i], "S_above", sep="_") ~ cld[row.names(cld)=="S_above",]$Letters,
                  conditions == paste(mod1$trait[i], "T_below", sep="_") ~ cld[row.names(cld)=="T_below",]$Letters,
                  conditions == paste(mod1$trait[i], "T_neutral", sep="_") ~ cld[row.names(cld)=="T_neutral",]$Letters,
                  conditions == paste(mod1$trait[i], "T_above", sep="_") ~ cld[row.names(cld)=="T_above",]$Letters,
                  conditions == paste(mod1$trait[i], "M_below", sep="_") ~ cld[row.names(cld)=="M_below",]$Letters,
                  conditions == paste(mod1$trait[i], "M_neutral", sep="_") ~ cld[row.names(cld)=="M_neutral",]$Letters,
                  conditions == paste(mod1$trait[i], "M_above", sep="_") ~ cld[row.names(cld)=="M_above",]$Letters))
              
              Tukey_results = all_traits_tables %>% drop_na(TukeyHSD)
              
              assign(as.vector(paste("Tukey",i,sep="_")), Tukey_results)
}

Tukey_results_all = rbind(Tukey_1,
                          Tukey_2,
                          Tukey_3,
                          Tukey_4,
                          Tukey_5,
                          Tukey_6,
                          Tukey_7,
                          Tukey_8,
                          Tukey_9,
                          Tukey_10,
                          Tukey_11,
                          Tukey_12,
                          Tukey_13,
                          Tukey_14)

# Create bubble chart #
library(ggplot2)
library(multcompView)
library(RColorBrewer)

level_order <- c('Above', 'Neutral', 'Below')
trait_order <- c('Aerobic', 'Anaerobic', 'Contains_Mobile_Elements', 'Facultatively_Anaerobic', 'Forms_Biofilms', 'Gram_Negative', 'Gram_Positive', 'Potentially_Pathogenic', 'Stress_Tolerant', 'M00175_Nitrogen_fixation_nitrogen_ammonia_, M00176_Sulfur_reduction_sulfate_H2S_', 'M00453_QseC_QseB_quorum_sensing_two_component_regulatory_system_', 'M00506_CheA_CheYBV_chemotaxis_two_component_regulatory_system_', 'M00513_LuxQN_CqsS_LuxU_LuxO_quorum_sensing_two_component_regulatory_system_')

Tukey_results_all$compartment <- factor(Tukey_results_all$compartment, levels=c('Mucus', 'Tissue', 'Skeleton'))

myPalette <- colorRampPalette(rev(brewer.pal(3, "Spectral")))

x = ggplot(Tukey_results_all, aes(x=factor(model_fit, level = level_order), trait, z = Z)) +
  #geom_point(aes(x = factor(model_fit, level = level_order))) +
  stat_summary_2d(aes(size = after_stat(value), color=after_stat(value)), geom = "point") +
  #scale_size_continuous(range = c(3, 10)) +
  labs(
  #  title = "Aerobic", 
     x = "Compartment & Model Fit", 
     y = "Functional Traits",
     size = "Z-score",
     color = "Z-score") +
    facet_wrap(~ compartment, strip.position="bottom") +
    theme_classic() +
  geom_text(data = Tukey_results_all, aes(x = model_fit, y = trait, label = TukeyHSD), size = 4, color = "black", hjust = -0.4, vjust = -0.8, fontface = "italic", check_overlap = TRUE) +
  scale_y_discrete(limits=rev) +
  scale_colour_gradient2(low="darkcyan", mid="grey", high="purple4")

x



------
  mod1 = all_comp_fit %>% group_by(compartment) %>% do(model = lm(Z ~ model_fit, data = .))

Letters1 <- multcompLetters4(mod1$model[[1]], TukeyHSD(aov(mod1$model[[1]])))
cld = as.data.frame.list(Letters1$model_fit)
all_comp_fit = all_comp_fit %>%
  mutate(TukeyHSD = case_when(
    compartment == "Mucus" & model_fit == "Below" ~ cld$Letters[1],
    compartment == "Mucus" & model_fit == "Neutral" ~ cld$Letters[2],
    compartment == "Mucus" & model_fit == "Above" ~ cld$Letters[3]))

Tukey_results = as.data.frame(all_comp_fit$sampleID)
colnames(Tukey_results) = 'sampleID'
Tukey_results$TukeyHSD_mucus = all_comp_fit$TukeyHSD

Letters2 <- multcompLetters4(mod1$model[[2]], TukeyHSD(aov(mod1$model[[2]])))
cld = as.data.frame.list(Letters2$model_fit)
all_comp_fit = all_comp_fit %>%
  mutate(TukeyHSD = case_when(
    compartment == "Skeleton" & model_fit == "Below" ~ cld$Letters[1],
    compartment == "Skeleton" & model_fit == "Neutral" ~ cld$Letters[2],
    compartment == "Skeleton" & model_fit == "Above" ~ cld$Letters[3]))

Tukey_results$TukeyHSD_skeleton = all_comp_fit$TukeyHSD

Letters3 <- multcompLetters4(mod1$model[[3]], TukeyHSD(aov(mod1$model[[3]])))
cld = as.data.frame.list(Letters3$model_fit)
all_comp_fit = all_comp_fit %>%
  mutate(TukeyHSD = case_when(
    compartment == "Tissue" & model_fit == "Below" ~ cld$Letters[1],
    compartment == "Tissue" & model_fit == "Neutral" ~ cld$Letters[2],
    compartment == "Tissue" & model_fit == "Above" ~ cld$Letters[3]))

Tukey_results$TukeyHSD_tissue = all_comp_fit$TukeyHSD

Tukey_results = Tukey_results %>% mutate(TukeyHSD = coalesce(TukeyHSD_mucus, TukeyHSD_skeleton, TukeyHSD_tissue)) %>% select(sampleID, TukeyHSD)

all_comp_fit$TukeyHSD = Tukey_results$TukeyHSD
