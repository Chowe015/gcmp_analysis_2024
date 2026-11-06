#load libraries
library(BiocManager)
library(dplyr)   
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(microbiome)
library(rstatix) # dunn_test
library(ANCOMBC)
library(ComplexHeatmap)# For advanced heatmaps
library(rngtools)
library(vegan)
library(RVAideMemoire)
library(multcompView)
library(purrr)
library(stringr)
library(patchwork)
library(readr)      # read_tsv(): fast reading of tab-separated files

## Step 2 Run Sloan Neutral Model to produce neutral table and figures
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
setwd("C:/Users/Colin Howe/Desktop/gcmp_final_analysis")

## Comparison between above neutral and below categories for mucus tissue and skeleton
## Import mapping and otu functional table for each compartment and create a specific phyloseq object
## Mucus functional phyloseq

sink_name = paste0("Picrust2_ancombc_analysis_log.txt")
sink(sink_name,append=FALSE,split=TRUE)

mucus_asv <- read.table(file = "./final_picrust2/mucus_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_mucus <- mucus_asv[, colSums(mucus_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
mucus_map <- read.table(file ="./final_picrust2/mucus_picrust2_out_pipeline/mucus_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
muc_otu <- otu_table(as.matrix(clean_mucus), taxa_are_rows = TRUE)

# Set metadata
muc_sample <- sample_data(mucus_map)

# Create Working phyloseq object
muc_funk_ph <- phyloseq(muc_otu,muc_sample)
muc_funk_ph

## Tissue functional phyloseq
tissue_asv <- read.table(file = "./final_picrust2/tissue_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_tissue <- tissue_asv[, colSums(tissue_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
tissue_map <- read.table(file ="./final_picrust2/tissue_picrust2_out_pipeline/tissue_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
tis_otu <- otu_table(as.matrix(clean_tissue), taxa_are_rows = TRUE)

# Set metadata
tis_sample <- sample_data(tissue_map)


# Create Working phyloseq object
tis_funk_ph <- phyloseq(tis_otu,tis_sample)
tis_funk_ph

## skeleton functional phyloseq
skeleton_asv <- read.table(file = "./final_picrust2/skeleton_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv",header=T,comment.char="", sep="\t", row.names = 1)
clean_skeleton <- skeleton_asv[, colSums(skeleton_asv != 0) != 0, drop = FALSE]

#### Import Metadata read.table
skeleton_map <- read.table(file ="./final_picrust2/skeleton_picrust2_out_pipeline/skeleton_combined_metadata.txt", header=T,comment.char="", sep="\t",row.names = 1)

# create matrix format for OTU and taxonomy table
ske_otu <- otu_table(as.matrix(clean_skeleton), taxa_are_rows = TRUE)

# Set metadata
ske_sample <- sample_data(skeleton_map)


# Create Working phyloseq object
ske_funk_ph <- phyloseq(ske_otu,ske_sample)
ske_funk_ph

## join each phyloseq object
merge_phyloseq(muc_funk_ph,tis_funk_ph,ske_funk_ph) -> merged_funk_ph

# Create data frame for downstream analysis
all_merge_df <- data.frame(sample_data(merged_funk_ph))

# join metacyc pathway descriptions as taxonomy 
taxonomy <- read.table(file = "./final_picrust2/metacyc_pathways_decrip_full.txt", sep = "\t", header = T,row.names = 1)
names(taxonomy)[1] <-paste("Phylum")
names(taxonomy)[2] <-paste("Class")
names(taxonomy)[3] <-paste("Order")

# Transform into a matrix
tax1 = tax_table(as.matrix(taxonomy))
tax_table(merged_funk_ph) <- tax1

merged_funk_ph

######## Raw Bubble plot Analysis ########################################
# Reorder the phyloseq object using prune_taxa.
## get relative abundance at specific functional level 
ps_abund <- merged_funk_ph %>% 
  tax_glom("Class") %>% 
  transform_sample_counts(function(x) x / sum(x))

## melt dataset with one taxa per row 
df <-psmelt(ps_abund)

# Non-parametric Test 
kw_test <- df %>% 
  group_by(Class,Phylum) %>%
  kruskal_test(Abundance ~ compartment_neutrality) %>% 
  adjust_pvalue(method="BH") %>% 
  add_significance("p.adj")

## subset only significance values 
sig_taxa <- kw_test %>% filter(p.adj < 0.05) %>% pull(Class)

# First Posthoc analysis 
dunn <- df %>% 
  filter(Class %in% sig_taxa) %>%
     group_by(Class) %>% 
  dunn_test(Abundance ~ compartment_neutrality, p.adjust.method = "holm")
dunn


## Write results for abundance comparative analysis 
## Output data tables
write.table(dunn, file="./final_picrust2/abundance_comp_neutrality_dunn_test.tsv",sep="\t", row.names = FALSE)
write.table(kw_test, file="./final_picrust2/abundance_comp_neutrality_kruskal_wallis_test.tsv",sep="\t", row.names = FALSE)

## Develop Bubble Plot #### 

# Compact letter display per taxon, from Dunn's adjusted p-values
cld <- dunn %>%
  mutate(comparison = paste(group1, group2, sep = "-")) %>%
  group_by(Class) %>%
  group_modify(~ {
    p <- setNames(.x$p.adj, .x$comparison)          # named vector: "A-B" = p
    l <- multcompLetters(p, threshold = 0.05)$Letters
    tibble(Group = names(l), Letter = unname(l))
  }) %>%
  ungroup()


plot_df <- df %>%
  filter(Class %in% sig_taxa) %>%
  group_by(Class, compartment_neutrality) %>%
  summarise(mean_ab = mean(Abundance), .groups = "drop") %>%
  group_by(Class) %>%
  mutate(z = as.numeric(scale(mean_ab))) %>%
  ungroup() #%>%
  
## Rename Column name to merge with cld & order taxa by overall abundance
colnames(plot_df)[2] <- "Group"

plot_df %>% left_join(cld, by = c("Class", "Group")) %>%
  mutate(Class = fct_reorder(Class, mean_ab, .fun = sum)) %>% right_join(kw_test, by="Class")-> plot_df

## Set each group in a specific order
plot_df$Group = factor(plot_df$Group, 
                              levels=c("mucus_above",
                                       "mucus_neutral",
                                       "mucus_below",
                                       "tissue_above",
                                       "tissue_neutral",
                                       "tissue_below",
                                       "skeleton_above",
                                       "skeleton_neutral",
                                       "skeleton_below"))
# 3. Bubble plot with significance letters
many_col <- c("#006CD1","black","#6A3D9A","#FF7F00","#E7298A","#E31A1C","#1B9E77","#E31A1C","#6A3D9A", "#FF7F00","#1B9E77","#A6CEE3")

cat_cols <- setNames(
  rep(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999"),
      length.out = n_distinct(plot_df$Category)),
  sort(unique(plot_df$Category)))

# map each label to its phylum color
lab_cols <- plot_df %>% distinct(Phylum,Class) %>% 
  mutate(col= as.character(Phylum)) 
lab_lookup <- setNames(lab_cols$Phylum, lab_cols$Class)


# 3. Bubble plot with significance letters
kw_bubble_plot <-ggplot(plot_df, aes(x = Group, y = Class)) +
  geom_point(aes(size = mean_ab * 100, fill = Phylum),
             shape = 21, colour = "grey30", stroke = 0.3) +
  geom_text(aes(label = Letter), size = 3.5, fontface = "bold",
            colour = "black", nudge_x = 0.28) +
  scale_size_continuous(range = c(0.5, 6), trans = "sqrt",
                        breaks = c(0.1, 1, 5, 20),
                        name = "Mean abundance (%)")  +
  scale_fill_manual(values = many_col, guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.7)) +
labs(x = NULL, y = NULL, caption = paste("Letters: Dunn's test (Holm-adjusted, p < 0.05) within each taxon.","Groups sharing a letter do not differ.")) +   
   facet_grid(Phylum ~ ., scales = "free_y", space = "free_y", switch ="y",
              labeller = labeller(Phylum = label_wrap_gen(18))) +
  scale_y_discrete(labels = lab_cols) + 
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10),
        plot.caption = element_text(hjust = 0))

kw_bubble_plot

ggsave("./final_picrust2/kw_dunn_bubble_plot.jpg", kw_bubble_plot, width = 20, height = 10, dpi= 600)

######################## Ancombc analysis ##############################################
# comparison between Neutral and non-neutral taxa
print(paste("Running ancombc2 analysis on Merged Phyloseq Object"))
set.seed(111)

# Set groups to as.factor
sample_data(merged_funk_ph)$compartment_neutrality <- as.factor(sample_data(merged_funk_ph)$compartment_neutrality)

#sample_data(merged_funk_ph)$compartment_neutrality <- relevel(sample_data(merged_funk_ph)$compartment_neutrality, ref="tissue_neutral")


## Run Ancombc2
output = ancombc2(data = merged_funk_ph, tax_level = "Class", fix_formula = "compartment_neutrality-1", group="compartment_neutrality", p_adj_method = "holm", global=TRUE, pairwise=TRUE,dunnet = TRUE)

## Pull out results
res_global = output$res_global
res = output$res
pair_res= output$res_pair

## Output data tables
write.table(res, file="./final_picrust2/ancombc2_picrust2_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)
write.table(pair_res, file="./final_picrust2/ancombc2_picrust2_pairwise_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)
write.table(res_global, file="./final_picrust2/ancombc2_picrust2_global_results_comp_neutrality.tsv",sep="\t", row.names = FALSE)


print(paste("Organizing results and generating Bubble Plot"))
df_fig_neutral = pair_res %>%
      dplyr::filter(diff_compartment_neutralitymucus_above == 1 |diff_compartment_neutralitymucus_below == 1 |diff_compartment_neutralitymucus_neutral == 1 |diff_compartment_neutralitytissue_above == 1 |diff_compartment_neutralitytissue_below == 1 |diff_compartment_neutralityskeleton_above == 1 |diff_compartment_neutralityskeleton_below == 1 |diff_compartment_neutralityskeleton_neutral == 1) %>% #diff_Intercept == 1 
  dplyr::mutate(lfc1 = ifelse(diff_compartment_neutralitymucus_above == 1, 
                              round(lfc_compartment_neutralitymucus_above, 2), 0),
                lfc2 = ifelse(diff_compartment_neutralitymucus_below == 1, 
                              round(lfc_compartment_neutralitymucus_below, 2), 0),
                lfc3 = ifelse(diff_compartment_neutralitymucus_neutral == 1, 
                              round(lfc_compartment_neutralitymucus_neutral, 2), 0),
                lfc4 = ifelse(diff_compartment_neutralitytissue_above == 1, 
                              round(lfc_compartment_neutralitytissue_above, 2), 0),
                lfc5 = ifelse(diff_compartment_neutralitymucus_below == 1, 
                              round(lfc_compartment_neutralitytissue_below, 2), 0),
                #lfc6 = ifelse(diff_compartment_neutralitytissue_neutral == 1, 
                 #             round(lfc_compartment_neutralitytissue_neutral, 2), 0),
                lfc7 = ifelse(diff_compartment_neutralityskeleton_above == 1, 
                              round(lfc_compartment_neutralityskeleton_above, 2), 0),
                lfc8 = ifelse(diff_compartment_neutralityskeleton_below == 1, 
                              round(lfc_compartment_neutralityskeleton_below, 2), 0),
                lfc9 = ifelse(diff_compartment_neutralityskeleton_neutral == 1, 
                              round(lfc_compartment_neutralityskeleton_neutral, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc9, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_neutral$group = recode(df_fig_neutral$group, 
                              `lfc1` = "mucus_above",
                              `lfc2` = "mucus_below",
                              `lfc3` = "mucus_neutral",
                              `lfc4` = "tissue_above",
                              `lfc5` = "tissue_below",
                              #`lfc6` = "tissue_neutral",
                              `lfc7` = "skeleton_above",
                              `lfc8` = "skeleton_below",
                              `lfc9` = "skeleton_neutral")

df_fig_neutral$group = factor(df_fig_neutral$group, 
                              levels=c("mucus_above",
                                       "mucus_neutral",
                                       "mucus_below",
                                       "tissue_above",
                                       #"tissue_neutral",
                                       "tissue_below",
                                       "skeleton_above",
                                       "skeleton_neutral",
                                       "skeleton_below"))



## Pairwise between neutrality and compartment
## Multicomp tool 
#taxon_name = "Glycan Biosynthesis"
## function that pulls q_values and applies multicomp tool analysis
get_sig_group <- function(taxon_name,pair_res){

## Filter data frame for each taxon
taxon_data <- pair_res %>% filter(taxon==taxon_name)

## Extract Column Names
q_vals <- taxon_data %>% 
  select(starts_with("q_"))

comp_names <-colnames(q_vals) %>% 
  gsub("q_","",.) %>%
  sub("_", "-", .)
  
p_vector <- as.numeric(q_vals[1,])
names(p_vector) <- comp_names

letters_obj <- multcompLetters(p_vector, threshold=0.05)
letters_obj

## Add significant label
#breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
#labels <- c("***", "**", "*", ".", "ns")

#sig_group <- cut(p_vector, breaks = breaks, label=labels, rights=TRUE)
#names(sig_group) <-names(p_vector)

# 5. Extract the letters array
letters_df <- data.frame(
  group = names(letters_obj$Letters),
  taxon=taxon_name,
  letter = letters_obj$Letters,  row.names = NULL) %>%
  dplyr::arrange(taxon_name)
return(letters_df)

#result_table <- data.frame(  group = names(p_vector),
#                             taxon=taxon_name,
#                             adjp = p_vector,
#                           sig_group= sig_group, row.names = NULL )%>%
#  dplyr::arrange(taxon_name)
#return(result_table)

#inner_join(result_table,letters_df, by= c("group","taxon"))
}

## apply to each pathway
all_taxa <- unique(pair_res$taxon)
all_letters_list <- lapply(all_taxa, function(x) get_sig_group(x, pair_res))
all_letters_df <- do.call(rbind, all_letters_list)


## Clean up group names to help merge letters with figure results
clean_letter_df <- all_letters_df %>%
  mutate(group = gsub("neutrality", "", group, ignore.case = TRUE),)%>%
  mutate(group = gsub("compartment_", "", group, ignore.case = TRUE),) 

sample_names <- tax_table(merged_funk_ph) %>% as.data.frame()
colnames(sample_names)[2] <- "taxon"

## Right join data frames
sig_fig_neutral <- df_fig_neutral %>% inner_join(clean_letter_df, by =c("group","taxon"), relationship = "many-to-many") 
test_fig <- sig_fig_neutral %>% right_join(sample_names, by ="taxon") %>% na.omit()

write.table(sig_fig_neutral, file="./final_picrust2/ancombc2_merged_letters_results.tsv", row.names=TRUE, sep = "\t")

#### Bubble Figure
sig_neutral_bubble_plot = test_fig %>%
  ggplot(aes(x = group, y = taxon)) + theme_minimal() +
  geom_vline(xintercept=0, linetype= "dashed", color = "grey50", linewidth= 0.5) +
  geom_point(aes(size=abs(value), fill= value >=0), shape=21, color="black", alpha=0.8) +
  scale_size_continuous(
    range = c(1, 10), 
    name = "Absolute Value"
  ) + scale_fill_manual(
    values = c("TRUE" = "#006CD1", "FALSE" = "#994F00"),
    labels = c("FALSE" = "Negative Lfc", "TRUE" = "Positive Lfc"),
    name = "Logfold Change") + 
  geom_text(aes(label = letter), 
            nudge_x = 0.18,       # Pushes letters slightly to the right of the bubble
            vjust = 0.5, 
            color = "black", 
            size = 4) +   
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y", switch ="y",
             labeller = labeller(Phylum = label_wrap_gen(18))) +
  scale_y_discrete(labels = lab_cols) + 
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10),
        plot.caption = element_text(hjust = 0))

sig_neutral_bubble_plot

ggsave("./final_picrust2/ancombc_logfold_change_comp_neutrality_bubbleplot.jpg", sig_neutral_bubble_plot, width = 20, height = 10, dpi= 600)


######## Heatmap #######################
# Heatmap Script developed with help from Claude Code, Sonnet 5.5 medium
## STEP 1: Pull pairwise results into a new tibble dataframe
heat_res <- as_tibble(pair_res)

# The name of your grouping variable as it appears inside the column names,
var_prefix <- "compartment_neutrality"

## Step 2: Reshape dataset to long format 

long_data <- heat_res %>%
  pivot_longer(
    cols = -taxon,
    names_to = c(".value", "contrast"),
    names_pattern = "^(lfc|se|W|p|q|diff)_(.*)$"
  )

## STEP 3: Define Pairwise comparisons and names for each pair
heat_pair <- long_data %>%
  filter(str_count(contrast, "compartment_neutrality") == 2) %>% # Filters contrast when there are two compartment_neutrality values.
  # split the comparison name on the variable name -> "", group1, group2
  separate(contrast, into = c("blank", "group1", "group2"),
           sep = "compartment_neutrality") %>% ## separated by compartment_neutrality
  select(-blank) %>%
  mutate(group1 = str_remove(group1, "_$")) %>%          # remove trailing "_"
  # split each group name, e.g. "mucus_above" -> compartment "mucus", class "above"
  separate(group1, into = c("comp1", "neu1"), sep = "_", remove = FALSE) %>%
  separate(group2, into = c("comp2", "neu2"), sep = "_", remove = FALSE)  

## factors for each comparison between neutrality and compartment
# The tables below say how to label a comparison, and "sign" = -1 means the
# stored direction is the reverse of the one we want, so the LFC is flipped.

within_key <- tribble(
  ~neu1,     ~neu2,     ~label, ~sign,
  "below",   "neutral", "B-N",   1,
  "neutral", "below",   "B-N",  -1,
  "above",   "neutral", "A-N",   1,
  "neutral", "above",   "A-N",  -1,
  "below",   "above",   "B-A",   1,
  "above",   "below",   "B-A",  -1
)

across_key <- tribble(
  ~comp1,     ~comp2,     ~label, ~sign,
  "mucus",    "tissue",   "M-T",   1,
  "tissue",   "mucus",    "M-T",  -1,
  "skeleton", "tissue",   "S-T",   1,
  "tissue",   "skeleton", "S-T",  -1,
  "skeleton", "mucus",    "S-M",   1,
  "mucus",    "skeleton", "S-M",  -1
)

## STEP 4: Build the two sets of contrasts
# Heatmap A: neutrality classes compared WITHIN a compartment
#   (same compartment, different neutrality class)
within <- heat_pair %>%
  filter(comp1 == comp2) %>%
  inner_join(within_key, by = c("neu1", "neu2")) %>%
  mutate(lfc = lfc * sign,        # flip sign where needed
         block = comp1,           # panel title = the compartment
         view = "within")

# Heatmap B: compartments compared WITHIN a neutrality class
#   (same neutrality class, different compartment)
across <- heat_pair %>%
  filter(neu1 == neu2) %>%
  inner_join(across_key, by = c("comp1", "comp2")) %>%
  mutate(lfc = lfc * sign,
         block = neu1,            # panel title = the neutrality class
         view = "across")

# Stack both sets into one table and keep only the columns we need
plot_df <- bind_rows(within, across) %>%
  select(taxon, view, block, label, lfc,se, q, diff)

# Quick check: how many tests are significant in each heatmap?
# (diff is TRUE where q < 0.05)
plot_df %>% count(view, diff)
plot_df %>% count(block, diff)

write.table(plot_df, file="./final_picrust2/heatmap_comp_neutral_results.tsv", sep="\t", row.names = FALSE)
## STEP 3: Define higher-level pathway classes
# A named list: each name is a class, each vector holds the pathway names that
# belong to it. Edit this list to move a pathway to a different class.
biosynthesis <- tibble(class = "Biosynthesis", taxon = c(
  "Amino Acid Biosynthesis", "Amine and Polyamine Biosynthesis",
  "Aromatic Compound Biosynthesis", "Carbohydrate Biosynthesis",
  "Cell Structure Biosynthesis",
  "Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis",
  "Fatty Acid and Lipid Biosynthesis", "Glycan Biosynthesis",
  "Nucleoside and Nucleotide Biosynthesis", "Secondary Metabolite Biosynthesis",
  "Storage Compound Biosynthesis",
  "superpathway of S-adenosyl-L-methionine biosynthesis"))

degradation <- tibble(class = "Degradation, utilization, assimilation", taxon = c(
  "1,5-anhydrofructose degradation", "Alcohol Degradation",
  "Amine and Polyamine Degradation", "Amino Acid Degradation",
  "Aromatic Compound Degradation", "C1 Compound Utilization and Assimilation",
  "Carbohydrate Degradation", "Carboxylate Degradation",
  "Chlorinated Compound Degradation",
  "Degradation/Utilization/Assimilation - Other", "ethylmalonyl-CoA pathway",
  "Fatty Acid and Lipid Degradation", "Inorganic Nutrient Metabolism",
  "Nucleoside and Nucleotide Degradation", "Polymeric Compound Degradation",
  "Secondary Metabolite Degradation",
  "superpathway of D-glucarate and D-galactarate degradation"))

metabolism <- tibble(class = "Precursor metabolites and energy", taxon = c(
  "Fermentation", "glyoxylate cycle", "Glycolysis", "Pentose Phosphate Pathways",
  "Respiration", "TCA cycle",
  "superpathway of glycolysis and the Entner-Doudoroff pathway"))

macromolecule <- tibble(class = "Macromolecule processing",
                        taxon = "Nucleic Acid Processing")


# Combine into one lookup table; factor() fixes the order of the classes
class_map <- bind_rows(biosynthesis, degradation, metabolism, macromolecule) %>%
  mutate(class = factor(class, levels = c(
    "Biosynthesis", "Degradation, utilization, assimilation",
    "Precursor metabolites and energy", "Macromolecule processing")))


# Attach the class to every row; stop with an error if a pathway has no class
plot_df <- plot_df %>% left_join(class_map, by = "taxon") 

# Row order: inside each class, sort pathways from highest to lowest average
# LFC across the neutrality contrasts (heatmap A).
taxon_order <- plot_df %>%
  filter(view == "within") %>%
  group_by(class, taxon) %>%
  summarise(mean_lfc = mean(lfc), .groups = "drop") %>%
  arrange(class, desc(mean_lfc)) %>%  pull(taxon)

# Turn the text columns into ordered factors so ggplot keeps our order.
# (ggplot puts the first level at the bottom, so rev() puts it at the top.)
plot_df <- plot_df %>%
  mutate(
    taxon_f = factor(taxon, levels = rev(taxon_order)),
    block   = factor(block, levels = c("mucus", "tissue", "skeleton",
                                       "above", "below","neutral")),
    label   = factor(label, levels = c( "A-N","B-N", "B-A",
                                       "M-T", "S-M","S-T" ))#,
    # cap extreme values so a few large LFCs do not wash out the colors
    #lfc_capped = pmax(pmin(lfc, lfc_cap), -lfc_cap)
  )

# Turn the list into a two-column lookup table: pathway name -> class.
# unlist() flattens the vectors; rep() repeats each class name once per pathway;
# factor(levels = ...) keeps the classes in the order written above.
class_map <- tibble(
  taxon = unlist(classes, use.names = FALSE),
  class = factor(rep(names(classes), lengths(classes)), levels = names(classes))
)

## STEP 6: Draft both heatmaps 
heat_color <- scale_fill_gradient2(
  low = "#006CD1", mid = "white", high = "#994F00", midpoint = 0,
  limits = c(-lfc_cap, lfc_cap),
  name = "log fold change (first group minus second)")

# dot colour: white on dark tiles, black on light tiles
heat_dots <- scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"))

# split the plot: rows = pathway class, columns = block (compartment or neutrality)
heat_facet <- facet_grid(class ~ block, scales = "free", space = "free",
                         switch = "y", labeller = labeller(class = label_wrap_gen(18)))


# Overall Plot Theme
heat_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        strip.placement = "outside",
        strip.text.x = element_text(face = "bold", size = 11),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.x = element_text(size = 8))

# Heatmap A: neutrality contrasts within each compartment
p_within <- ggplot(filter(plot_df, view == "within"),
                   aes(x = label, y = taxon_f, fill = lfc)) +
  geom_tile(color = "white", linewidth = 0.3) +             # one tile per test
  geom_point(data = filter(plot_df, view == "within", diff),  # dot if q < 0.05
             aes(color = abs(lfc) > 1.8), size = 1.1, show.legend = FALSE) +
  heat_color + heat_dots + heat_facet +
  labs(title = "Neutrality contrasts within each compartment", x = NULL, y = NULL) +
  heat_theme+
theme(axis.text.y = element_blank(), strip.text.y.left = element_blank())


# Heatmap B: compartment contrasts within each neutrality class
# (same code, but pathway names and class labels are hidden because
#  they already appear on heatmap A)
p_across <- ggplot(filter(plot_df, view == "across"),
                   aes(x = label, y = taxon_f, fill = lfc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = filter(plot_df, view == "across", diff),
             aes(color = abs(lfc) > 1.8), size = 1.1, show.legend = FALSE) +
  heat_color + heat_dots + heat_facet +
  labs(title = "Compartment contrasts within each neutrality class", x = NULL, y = NULL) +
  heat_theme 
## STEP 7: Combine the two heatmaps and save
# "+" puts the plots side by side; guides = "collect" shares one color legend
final <- ( p_across +p_within) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = paste(
    "compartments: M = mucus, T = tissue,S = skeleton.",
     "neutrality: A = above, B = below, N = neutral, Dots: adjp < 0.05, holm.",
    "\n Each column = first group minus second.",
    "Warm colors = enrichment in first Group")) &
  theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"))

# width and height are in inches; use .pdf or .svg in outfile for vector output
ggsave("./final_picrust2/ancombc2_heatmap_lfc_plot.pdf", final, width = 15, height = 12, dpi = 600, bg = "white")

