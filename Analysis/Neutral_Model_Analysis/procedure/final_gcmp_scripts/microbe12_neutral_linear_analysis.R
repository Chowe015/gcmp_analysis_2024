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

sink_name = paste0("./alpha_beta_diversity/microbe12_gcmp_neutral_linear_model_log.txt")
sink(sink_name,append=FALSE,split=TRUE)


##  Compartment Boxplot with R2 
coral_color =c("#009197","#EC8C5C","#8E3B97")
miseq_table <- read.table(file = "./miseq_neutral_results_table.txt",header=T,sep="\t")
hiseq_table <- read.table(file = "./hiseq_neutral_results_table.txt",header=T,sep="\t")

# remove low count samples
filter(miseq_table, Samples>10 & Samples <60) -> miseq_table

print(paste("Miseq Datatable in Use"))
compartment = miseq_table$compartment
region = miseq_table$region
Rsquared = miseq_table$Rsqr
faith = miseq_table$faithpd
samples = miseq_table$Samples
region = factor(miseq_table$region)
expedition = miseq_table$expedition
migration = miseq_table$m

print(paste("Running linear regression analysis for Rsquared and Migration"))
rsqr_test <- lm(Rsquared ~ faith+compartment+region+expedition)
summary(rsqr_test) 

faith_rsqr_test <- lm(Rsquared ~ faith)
summary(faith_rsqr_test )

compartment_rsqr_test <- lm(Rsquared ~ compartment)
summary(compartment_rsqr_test)

expedition_rsqr_test <- lm(Rsquared ~ expedition)
summary(expedition_rsqr_test)

samples_rsqr_test <- lm(Rsquared ~ samples)
summary(samples_rsqr_test)

####################### migration comparisons ##########################################
## Migration 
migration_test <- lm(migration ~ faith+compartment+region+expedition)
summary(migration_test) 

faith_migration_test <- lm(migration ~ faith)
summary(faith_migration_test) 

compartment_migration_test <- lm(migration ~ compartment)
summary(compartment_migration_test)

expedition_migration_test <- lm(migration ~ expedition)
summary(expedition_migration_test)

####################### Rsquared comparisons ##########################################
print(paste("Plotting Linear Models and Boxplots for miseq_table dataset"))

## Plot Linear regression Host Faith Pd
faith_r_summary <-summary(faith_rsqr_test)
faith_r_pvalue <- paste0("pvalue: ",round(faith_r_summary$coefficients[2,4],3))

faith_r2_plot = plot(faith,Rsquared, pch=16,cex=1.3,col="#009197", main = "Rsquared and Coral Faith Pd", xlab = "Coral Species Diversity (faith Pd)", ylab ="Rquared (R2)") +
  text(x=1500, y=.7, label = faith_r_pvalue) + abline(lm(Rsquared ~ faith), col="#EC8C5C") 

## Compartment Boxplot
testp <- ggplot(miseq_table, aes(y=Rsquared, x=compartment)) + geom_boxplot(aes(color=compartment)) +
  theme_classic() + geom_point(aes(color=compartment)) +labs(xlab="R2", ylab="Coral Compartments", title ="Compartment Model Fit") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#009197","#EC8C5C","#8E3B97")) + geom_signif(comparisons = list(c("mucus", "tissue"), c("mucus", "skeleton"),c("tissue", "skeleton")), y_position = c(.85,1, .75))
testp
## Plot Linear regression Samples
## Retrieve and assign the pvalue from summary
samp_r_summary <-summary(samples_rsqr_test)
samp_r_pvalue <- paste0("pvalue: ",round(samp_r_summary$coefficients[2,4],3))

samples_r2_plot = plot(samples,Rsquared, pch=16,cex=1.3,col="#009197", main = "Rsquared by Expedition by Compartment Sample Size", xlab = "Compartment/Expedition Sample Count (n)", ylab ="Rquared (R2)")+ text(x=45, y=0.2, label = samp_r_pvalue) + abline(lm(Rsquared ~ samples), col="#EC8C5C") 

## Plot Linear regression migration 
###################### Rsquared comparisons ##########################################
## Plot Linear regression Host Faith Pd
faith_m_summary <-summary(faith_migration_test)
faith_m_pvalue <- paste0("pvalue: ",round(faith_m_summary$coefficients[2,4],3))

faith_m_plot = plot(faith,migration, pch=16,cex=1.3,col="#009197", main = "Migration and Coral Faith Pd", xlab = "Coral Species Diversity (faith Pd)", ylab ="Migration (m)") + text(x=2000, y=.07, label = faith_m_pvalue) + abline(lm(migration ~ faith), col="#EC8C5C") 

## Plot Linear regression Compartment
## Assign levels 
compartment <- factor(compartment, levels = c("mucus", "tissue", "skeleton")) 
# Compartment Boxplot
compart_m_plot = plot(compartment,migration,col=coral_color, main = "Migration and Coral Anatomical Comartments", xlab = "Coral Compartments", ylab ="Migration (m)")+
  scale_color_manual(values=c("#009197","#EC8C5C","#8E3B97")) + geom_signif(comparisons = list(c("Mucus", "Tissue"), c("Mucus", "Skeleton"),c("Tissue", "Skeleton")), y_position = c(.85,1, .75))
compart_m_plot
## Plot Linear regression Samples
## Retrieve and assign the pvalue from summary
samp_s_summary <-summary(expedition_migration_test)
samp_s_pvalue <- paste0("pvalue: ",round(samp_s_summary$coefficients[2,4],3))

samples_m_plot = plot(samples,migration, pch=16,cex=1.3,col="#009197", main = "Migration by Sample Size", xlab = "Compartment/Expedition Sample Count (n)", ylab ="Migration (m)") + text(x=45, y=.06, label = samp_s_pvalue) + abline(lm(migration~samples), col="#EC8C5C") 

## Output Alpha Diversity Plot plots as Pdf.
print(paste("Exporting linear Regression and Box Plots as pdf"))

# Migration & Sample Size
sample_m_name <- paste0("./alpha_beta_diversity/miseq_table_migration_sample_linear_model.pdf")
#ggsave(filename = sample_m_name, plot=samples_m_plot)

# Migration & Coral Faith Pd
faith_m_name <- paste0("./alpha_beta_diversity/miseq_table_migration_faith_linear_model.pdf")
#ggsave(filename=faith_m_name,faith_m_plot)#width = 20, height = 10, dpi= 300)

# Migration & Coral Compartment
comp_m_name <- paste0("./alpha_beta_diversity/miseq_table_migration_compartment_boxplot.pdf")
#ggsave(compart_m_plot, filename=comp_m_name)

# Rsquared & Sample Size
sample_r_name <- paste0("./alpha_beta_diversity/miseq_table_Rsquared_sample_linear_model.pdf")
#ggsave(samples_r2_plot, filename=sample_r_name)

# Rsquared & Coral Faith Pd
faith_r_name <- paste0("./alpha_beta_diversity/miseq_table_Rsquared_faith_linear_model.pdf")
#ggsave(faith_r2_plot, filename=faith_r_name)

# Rsquared & Coral Compartment
comp_r_name <- paste0("./alpha_beta_diversity/miseq_table_Rsquared_compartment_boxplot.pdf")
#ggsave(testp, filename=comp_r_name)

####### Hiseq ########
print(paste("Hiseq Datatable in Use"))
compartment = hiseq_table$compartment
region = hiseq_table$region
Rsquared = hiseq_table$Rsqr
faith = hiseq_table$faithpd
samples = hiseq_table$Samples
region = factor(hiseq_table$region)
expedition = hiseq_table$expedition
migration = hiseq_table$m

print(paste("Running linear regression analysis for Rsquared and Migration"))
rsqr_test <- lm(Rsquared ~ faith+compartment+region+expedition)
summary(rsqr_test) 

faith_rsqr_test <- lm(Rsquared ~ faith)
summary(faith_rsqr_test )

compartment_rsqr_test <- lm(Rsquared ~ compartment)
summary(compartment_rsqr_test)

expedition_rsqr_test <- lm(Rsquared ~ expedition)
summary(expedition_rsqr_test)

samples_rsqr_test <- lm(Rsquared ~ samples)
summary(samples_rsqr_test)

####################### migration comparisons ##########################################
## Migration 
migration_test <- lm(migration ~ faith+compartment+region+expedition)
summary(migration_test) 

faith_migration_test <- lm(migration ~ faith)
summary(faith_migration_test) 

compartment_migration_test <- lm(migration ~ compartment)
summary(compartment_migration_test)

expedition_migration_test <- lm(migration ~ expedition)
summary(expedition_migration_test)

####################### Rsquared comparisons ##########################################
print(paste("Plotting Linear Models and Boxplots for hiseq_table dataset"))

## Plot Linear regression Host Faith Pd
faith_r_summary <-summary(faith_rsqr_test)
faith_r_pvalue <- paste0("pvalue: ",round(faith_r_summary$coefficients[2,4],3))

faith_r2_plot = plot(faith,Rsquared, pch=16,cex=1.3,col="#009197", main = "Rsquared and Coral Faith Pd", xlab = "Coral Species Diversity (faith Pd)", ylab ="Rquared (R2)") + 
  text(x=1500, y=.7, label = faith_r_pvalue) + abline(lm(Rsquared ~ faith), col="#EC8C5C") 

## Compartment Boxplot
testp <- ggplot(hiseq_table, aes(y=Rsquared, x=compartment)) + geom_boxplot(aes(color=compartment)) +
  theme_classic() + geom_point(aes(color=compartment)) +labs(xlab="R2", ylab="Coral Compartments", title ="Compartment Model Fit") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#009197","#EC8C5C","#8E3B97")) + geom_signif(comparisons = list(c("Mucus", "Tissue"), c("Mucus", "Skeleton"),c("Tissue", "Skeleton")), y_position = c(.85,1, .75))

## Plot Linear regression Samples
## Retrieve and assign the pvalue from summary
samp_r_summary <-summary(samples_rsqr_test)
samp_r_pvalue <- paste0("pvalue: ",round(samp_r_summary$coefficients[2,4],3))

samples_r2_plot = plot(samples,Rsquared, pch=16,cex=1.3,col="#009197", main = "Rsquared by Expedition by Compartment Sample Size", xlab = "Compartment/Expedition Sample Count (n)", ylab ="Rquared (R2)")+ text(x=50, y=0.7, label = samp_r_pvalue) + abline(lm(Rsquared ~ samples), col="#EC8C5C") 

## Plot Linear regression migration 
###################### Rsquared comparisons ##########################################
## Plot Linear regression Host Faith Pd
faith_m_summary <-summary(faith_migration_test)
faith_m_pvalue <- paste0("pvalue: ",round(faith_m_summary$coefficients[2,4],3))

faith_m_plot = plot(faith,migration, pch=16,cex=1.3,col="#009197", main = "Migration and Coral Faith Pd", xlab = "Coral Species Diversity (faith Pd)", ylab ="Migration (m)") + text(x=1500, y=.07, label = faith_m_pvalue) + abline(lm(migration ~ faith), col="#EC8C5C") 

## Plot Linear regression Compartment
## Assign levels 
compartment <- factor(compartment, levels = c("mucus", "tissue", "skeleton")) 
# Compartment Boxplot
compart_m_plot = plot(compartment,migration,col=coral_color, main = "Migration and Coral Anatomical Comartments", xlab = "Coral Compartments", ylab ="Migration (m)") 

## Plot Linear regression Samples
## Retrieve and assign the pvalue from summary
samp_s_summary <-summary(expedition_migration_test)
samp_s_pvalue <- paste0("pvalue: ",round(samp_s_summary$coefficients[2,4],3))

samples_m_plot = plot(samples,migration, pch=16,cex=1.3,col="#009197", main = "Migration Sample Size", xlab = "Compartment/Expedition Sample Count (n)", ylab ="Migration (m)") + text(x=55, y=.07, label = samp_s_pvalue) + abline(lm(migration~samples), col="#EC8C5C") 

## Output Alpha Diversity Plot plots as Pdf.
print(paste("Exporting linear Regression and Box Plots as pdf"))

# Migration & Sample Size
#sample_m_name <- paste0("./alpha_beta_diversity/hiseq_table_migration_sample_linear_model.pdf")
#ggsave(samples_m_plot, filename=sample_m_name)

# Migration & Coral Faith Pd
#faith_m_name <- paste0("./alpha_beta_diversity/hiseq_table_migration_faith_linear_model.pdf")
#ggsave(faith_m_plot, filename=faith_m_name)

# Migration & Coral Compartment
#comp_m_name <- paste0("./alpha_beta_diversity/hiseq_table_migration_compartment_boxplot.pdf")
#ggsave(compart_m_plot, filename=comp_m_name)

# Rsquared & Sample Size
#sample_r_name <- paste0("./alpha_beta_diversity/hiseq_table_Rsquared_sample_linear_model.pdf")
#ggsave(samples_r2_plot, filename=sample_r_name)

# Rsquared & Coral Faith Pd
#faith_r_name <- paste0("./alpha_beta_diversity/hiseq_table_Rsquared_faith_linear_model.pdf")
#ggsave(faith_r2_plot, filename=faith_r_name)

# Rsquared & Coral Compartment
#comp_r_name <- paste0("./alpha_beta_diversity/hiseq_table_Rsquared_compartment_boxplot.pdf")
#ggsave(testp, filename=comp_r_name)

print(paste("Finished!"))

