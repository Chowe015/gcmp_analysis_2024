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


#Import from .qza file into a phyloseq object
asv <- qza_to_phyloseq(features = "./hiseq_files_2018/feature_table_decon_all_1000.qza")

#### Import Metadata read.table
metadata <- read.table(file = "./hiseq_files_2018/GCMP_EMP_map_r29.txt",header=T,comment.char="", row.names=1, sep="\t")

### Import Tree file from biom output tree.nwk
tree <- read_tree("./hiseq_files_2018/tree.nwk")

### Import taxonomy from biom output as .tsv format using read.table
taxonomy <- read.table(file = "./hiseq_files_2018/taxonomy.tsv", sep = "\t", header = T ,row.names = 1)


# clean the taxonomy
##code referenced from Yan Hui: email me@yanh.org github: yanhui09
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

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

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

###Agglomerate taxa to family 
glom <- tax_glom(subject, taxrank = 'Family', NArm = TRUE)
glom

phyloseq::tax_table(glom)%>%
  as.data.frame()%>%
  rownames_to_column("id") -> glom_tax

#Bruns et al. 2016 Neutral Model functions
################define sncm.fit############
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p),method="Nelder-Mead")
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.01, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p_abundance', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}
##Identify above, below and neutral OTUs 

addf = function(freq, pred.upr, pred.lwr, padj) {
  if(freq > pred.upr){
    model = "above"
  }
  else if(freq < pred.lwr){
    model = "below"
  }
  else if(pred.upr >= freq & freq >= pred.lwr){
    model = "neutral"
  }
  return(model)
}


table <- c("M","T","S")
for (i in table){
print(paste("Run Neutral Model on",i," OTU Table"))

subset_samples(glom,tissue_compartment == i) -> subject

  
paste(print("Creating glom mapping file"))
phyloseq::sample_data(subject) %>% as.data.frame() -> glom_map
    
biosample <- glom_map$sample_type[2]
  
###coral communities
set.seed(111)
spp <- prune_taxa(taxa_sums(subject) > 0, subject)
spp
otu_spp = as.data.frame(t(otu_table(spp)))

neutral_mod_w = sncm.fit(otu_spp)
neutral_mod1_w = sncm.fit(otu_spp,stats = F)
print(length(neutral_mod1_w$freq))

neutral_mod1_w$model = mapply(addf, neutral_mod1_w$freq, neutral_mod1_w$pred.upr, neutral_mod1_w$pred.lwr)  
print(length(which(neutral_mod1_w$model == "above")))
print(length(which(neutral_mod1_w$model == "below")))
print(length(which(neutral_mod1_w$model == "neutral")))

neutral_mod1_w %>% mutate(ci.mean = (pred.upr-pred.lwr),
                          SE = (ci.mean/(2*1.96)),
                          zstat = (freq/SE),
                          p_value = (exp(-0.717*zstat-0.416*zstat^2)),
                          padj = (exp(-0.717*zstat-0.416*zstat^2))*length(neutral_mod1_w$model)) %>%
  as.data.frame() %>% rownames_to_column("id")-> neutral_results_padj

print(paste("joining neutral table and ref taxonomy"))
## Inner_join feature table with neutral model families
inner_join(glom_tax,neutral_results_padj, by = "id") %>% as.data.frame() -> nonneutral_tax

# Print output files
print(paste("Writing Neutral Model tsv"))
neutral_table_name <- paste0(biosample,"_neutral_model.tsv")
write.table(nonneutral_tax,neutral_table_name,row.names = FALSE,sep="\t", col.names = TRUE)

#################### plot neuModel ################################
r <- round(neutral_mod_w$Rsqr,3)
m_var <- round(neutral_mod_w$m,3)
lb1 <- paste0("Rsquared: ",r)
lb2 <- paste0("migration (m) :",m_var)
lb3 <- paste0(biosample," Neutral_Model")
lb4 <-paste0(biosample," Log10 (Mean Relative abundance)")
lb5 <- paste0(biosample," Occurance frequency")

###### draft neutral model plot #############
p_neuM_coral=ggplot(neutral_results_padj, aes(x=log10(p_abundance), y=freq, color = model))+
  geom_point(shape=19, alpha=0.7, size=2.5)+
  geom_line(aes(x=log10(p_abundance), y=freq.pred),
            color="blue",linetype="solid",linewidth=1)+
  geom_line(aes(x=log10(p_abundance), y=pred.lwr),
            color="blue",linetype="dashed",linewidth=1)+
  geom_line(aes(x=log10(p_abundance), y=pred.upr),
            color="blue",linetype="dashed",linewidth=1)+
  scale_colour_manual(name="OTUs",
                      values=c("#E41A1C", "#999999", "#FF7F00","purple"),
                      breaks=c("above", "neutral", "below","neutral_significant"),
                      labels=c("Over-represented", "Neutrally-distributed", "Under-represented","Significant_low_weight"))+
  xlab(lb4) + ylab(lb5)+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=12), 
        legend.title = element_text(size=11.5, face = "bold"),
        legend.text = element_text(size=11.5),
        legend.position="none") + ggtitle(lb3)+ theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"))

p_neuM_coral1 = p_neuM_coral + annotate("text", x = -4.0, y=0.9, size=6, label=lb1, 
                                        fontface="bold", color="black",
                                        parse=TRUE) +
  annotate("text", x = -4.0, y=0.8, size=6, label=lb2, 
           fontface="bold", color="black",
           parse=TRUE)


print(paste("Creating Neutral Figure pdf files"))
## Print Neutral Model Plot
neutral_file_name <- paste0(biosample,"_neutral_model_Plot.pdf")
ggsave(p_neuM_coral1, filename=neutral_file_name,width = 20, height = 10, dpi= 600)

#### create ASV tables by id ** This file will be used in microbe_neutral_compartment.R and picrust2_neutral_table_generator.R
print(paste("Generating Agglomerated ASV Table dataset..."))
phyloseq::otu_table(subject)%>%
  as.data.frame()%>%
  rownames_to_column("id") -> glom_otu_table

## Output .tsv from the otu table file
glom_otu_name <- paste0(i,"_glom_table.tsv")
write.table(glom_otu_table, file =glom_otu_name ,sep = "\t",row.names = FALSE)

#### create taxonomy tables by id  ** This 
paste(print("Printing Agglomerated Taxonomy Table"))

phyloseq::tax_table(subject)%>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  select(-c("Genus","Species"))-> glom_taxonomy

## Output .tsv from the taxonomy file
glom_taxonomy_name <- paste0(i,"_glom_taxonomy.tsv")
write.table(glom_taxonomy, file =glom_taxonomy_name,sep = "\t", row.names = FALSE)
## Output .csv from the taxonomy file
glom_mapping_name <- paste0(i,"_glom_metadata.tsv")
write.table(glom_map, file =glom_mapping_name,sep = "\t", row.names = TRUE)

}

print(paste("Finished!"))
