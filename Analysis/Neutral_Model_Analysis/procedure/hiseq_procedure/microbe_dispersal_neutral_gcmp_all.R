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
library(bbmle)
library(biomformat)


#Step 1: Import qiime2 tables, mapping, taxonomy, tree and coral phylogeny across mucus tissue and skeleton
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)
feature_table_path <-args[1]
metadata_path <-args[2]
taxonomy_path <-args[3]
tree_path <-args[4]


sink_name = paste0("step_5_GCMP_Miseq_Dispersal_compartments_log.txt")
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

#### Import Tree file from biom output tree.nwk

#print(paste("tree_path",tree_path))
tree <- read_tree(tree_path)



#### Import taxonomy from biom output as .tsv format using read.table

print(paste("Loading Taxonomy text files from path:", taxonomy_path))
taxonomy <- read.table(file = taxonomy_path, sep = "\t", header = T ,row.names = 1)

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

# Inital Subset for agglomeration 
sub <- subset_samples(phylo, outgroup=="n" & subset_col=="yes")
sub

### rarefy to even depth
print(paste("Generating Rarefied Coral dataset..."))
rarefied = rarefy_even_depth(sub, rngseed=111, sample.size=1000, replace=F, trimOTUs = TRUE)
print(rarefied)

#### Agglomerate taxa to family 
print(paste("Agglomerated Taxonomy to the Family Level"))
glom <- tax_glom(rarefied, taxrank = 'Family', NArm = TRUE)
print(glom)

print(paste("Create Neutral Model Function"))
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


print(paste("Starting Loop across expeditions"))
################# Start for loop across expeditions ######################
##notes there are no skeleton or tissue samples from Hiseq E8 dataset
group_a <- c("M","S","T")
group_b <- c("S","T","M")
for (i in group_a){
  for (g in group_b){
# Subset only corals from database
spp_subset <-subset_samples(main_phylo, tissue_compartment == print(paste(i)))
pool_subset <-subset_samples(main_phylo, tissue_compartment == print(paste(g)))

spp_df <- sample_data((spp_subset))
spp_name = spp_df$BiologicalMatter[2]

pool_df <- sample_data((pool_subset))
pool_name = pool_df$BiologicalMatter[2]

###coral communities
set.seed(111)
spp <- prune_taxa(taxa_sums(spp_subset) > 0, spp_subset)
spp
otu_spp = as.data.frame(t(otu_table(spp)))

pool <- prune_taxa(taxa_sums(pool_subset) > 0, pool_subset)
pool

otu_pool = as.data.frame(t(otu_table(pool)))

## Neutral fit
mod_w_coral = sncm.fit(otu_spp,otu_pool)
mod1_w_coral = sncm.fit(otu_spp,otu_pool,stats = F)

mod1_w_coral %>% mutate(ci.mean = (pred.upr-pred.lwr),
                        SE = (ci.mean/(2*1.96)),
                        zstat = (freq/SE),
                        p_value = (exp(-0.717*zstat-0.416*zstat^2)),
                        padj = (exp(-0.717*zstat-0.416*zstat^2))*length(mod1_w_coral$model)) %>%
  as.data.frame() %>% rownames_to_column("id")-> neutral_results_padj

neutral_results_padj$model = mapply(addf, neutral_results_padj$freq, neutral_results_padj$pred.upr, neutral_results_padj$pred.lwr, neutral_results_padj$padj)  

length(which(neutral_results_padj$model == "above"))
length(which(neutral_results_padj$model == "below"))
length(which(neutral_results_padj$model == "neutral"))


neutral_results_padj$model <- factor(neutral_results_padj$model, 
                                     levels = c("above", "neutral", "below"))

## Inner_join feature table with neutral model families
inner_join(neutral_results_padj,glom_tax, by = "id") %>% as.data.frame() -> nonneutral_tax

psuedo_spp_name <- paste0(spp_subset,"_subset_",pool_subset,"_abundance_neutral_table.tsv")
(write.table(nonneutral_tax, file =psuedo_spp_name,sep="\t", row.names = FALSE))

############## Switch Frequency and Abundance######################################################
## Neutral fit
mod_w_coral1 = sncm.fit(otu_pool,otu_spp)
mod1_w_coral1 = sncm.fit(otu_pool,otu_spp,stats = F)

mod1_w_coral1 %>% mutate(ci.mean = (pred.upr-pred.lwr),
                         SE = (ci.mean/(2*1.96)),
                         zstat = (freq/SE),
                         p_value = (exp(-0.717*zstat-0.416*zstat^2)),
                         padj = (exp(-0.717*zstat-0.416*zstat^2))*length(mod1_w_coral1$model)) %>%
  as.data.frame() %>% rownames_to_column("id")-> neutral_results_padj2

neutral_results_padj2$model = mapply(addf, neutral_results_padj2$freq, neutral_results_padj2$pred.upr, neutral_results_padj2$pred.lwr, neutral_results_padj2$padj)  

length(which(neutral_results_padj2$model == "above"))
length(which(neutral_results_padj2$model == "below"))
length(which(neutral_results_padj2$model == "neutral"))


neutral_results_padj2$model <- factor(neutral_results_padj2$model, 
                                      levels = c("above", "neutral", "below"))

inner_join(neutral_results_padj2,glom_tax, by = "id") %>% as.data.frame() -> nonneutral_tax2
psuedo_pool_name <- paste0(pool_subset,"_subset_",spp_subset,"_abundance_neutral_table.tsv")
  (write.table(nonneutral_tax2, file =psuedo_pool_name,sep="\t", row.names = FALSE))
###Set Results parameters for plot statistics
r <- round(mod_w_coral$Rsqr,3)
m_var <- round(mod_w_coral$m,3)
lb1 <- paste0("R_squared",":",r)
lb2 <- paste0("migration_rate (m):",m_var)
lb3 <- paste0(spp_name,"&",pool_name," Disperal Neutral Model")
lb4 <-paste0(pool_name," Log10 (Mean Relative abundance)")
lb5 <- paste0(spp_name," Occurance frequency")

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

########################## Switch Frequency and Abundance Plot
#plot neuModel

r <- round(mod_w_coral1$Rsqr,3)
m_var <- round(mod_w_coral1$m,3)
lb1 <- paste0("R_squared",":",r)
lb2 <- paste0("migration_rate (m):",m_var)
lb3 <- paste0(paste0(spp_name,"&",pool_pool," Disperal Neutral Model"))
lb4 <-paste0(spp_name," Log10 (Mean Relative abundance)")
lb5 <- paste0(pool_name," Occurance frequency")

p_neuM_coral2=ggplot(neutral_results_padj2, aes(x=log10(p_abundance), y=freq, color = model))+
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


neutral_plot2 = p_neuM_coral2 + annotate("text", x = -4.0, y=0.9, size=6, label=lb1, 
                                         fontface="bold", color="black",
                                         parse=TRUE) +
  annotate("text", x =-4.0, y=0.8, size=6, label=lb2, 
           fontface="bold", color="black",
           parse=TRUE)

}
 }


print("Finished!")