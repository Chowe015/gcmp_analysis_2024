#Loading Required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(minpack.lm)
library(Hmisc)
library(bbmle)
library(biomformat)


## Step 2 Run Sloan Neutral Model to produce neutral table and figures
#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

glom_table_path <-args[1]
taxonomy_path <- args[2]
metadata_path <-args[3]
biosample <- args[4]

sink_name = paste0(biosample," ","Neutral_model_results_log.txt")
sink(sink_name,append=FALSE,split=TRUE)


#print(paste("Imported GCMP OTU Table"))
## import table 
glom_table <- read.table(glom_table_path, sep = "\t",header = TRUE,row.names=1,check.name=FALSE)
#import taxonomy
glom_tax <- read.table(taxonomy_path, sep = "\t",header = TRUE,check.name=FALSE)
# import mapping file
glom_mapping <-read.table(metadata_path, sep = ",",header = TRUE,check.name=FALSE)

## Testing different import tools
#glom_table <- read.table("./Tissue/T_glom_table.tsv", sep = "\t",header = TRUE,row.names=1,check.name=FALSE)
##import taxonomy
#glom_tax <- read.table("./Tissue/T_glom_taxonomy.tsv", sep = "\t",header = TRUE,check.name=FALSE)
# import mapping file
#glom_mapping <-read.table("./Tissue/T_glom_metadata.csv", sep = ",",header = TRUE,check.name=FALSE)


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
  
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p),"Nelder-Mead")
  
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

print(paste("Run Neutral Model Function on OTU Table"))
## pull in otu table 

otu_frame = as.data.frame(t(glom_table))
## Neutral fit
neutral_mod_w = sncm.fit(otu_frame)
neutral_mod1_w = sncm.fit(otu_frame,stats = F)


print(length(neutral_mod_w$freq))

neutral_mod1_w$model = mapply(addf, neutral_mod1_w$freq, neutral_mod1_w$pred.upr, neutral_mod1_w$pred.lwr)  
length(which(neutral_mod1_w$model == "above"))
length(which(neutral_mod1_w$model == "below"))
length(which(neutral_mod1_w$model == "neutral"))

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

############################### create compartment pseudo table outputs #########################
print(paste("Writing Above and Below Pseudo tables"))
glom_table2 <-read.table(glom_table_path, sep = "\t",header = TRUE,check.name=FALSE)
#glom_table2 <-read.table("./Tissue/T_glom_table.tsv", sep = "\t",header = TRUE,check.name=FALSE)
# Subset Taxonomy and join with neutral model results
neutral = c("above","below")
for (b in neutral){
  ## Inner_join will connect  ASV id with taxonomy using id column
  inner_join(neutral_results_padj,glom_tax, by = "id") -> non_neutral_tax
  
  ## Subset ASV from each model and significance 
  sig_neutral_tax = select(filter(non_neutral_tax, model == b & padj <=0.05),id,p_abundance,freq,padj,Phylum,Class,Order,Family,model)
  
  #subset otu table based on neutral results and rename column names
  neutral_table <-subset(glom_table2, id  %in%  sig_neutral_tax$id)
  l = length(neutral_table)
  colnames(neutral_table)[2:l] <- paste(colnames(neutral_table)[2:l], b, sep = "_")
 colnames(neutral_table)[1] <- "id" 
  ## Create new columns for the metadata data
  glom_mapping %>% mutate(psuedo_sample_name = paste(colnames(neutral_table)[2:l]),
                          compartment_neutrality = paste(biosample,b, sep = '_'),
                          compartment = paste (biosample),
                          neutrality = paste(b)) %>% 
    select(psuedo_sample_name,sample_name_backup, host_scientific_name,expedition_number,political_area,ocean_area,compartment_neutrality,compartment,neutrality) %>% 
    as.data.frame() -> glom_metadata
  
  ## Subset significant results with taxonomy to retain only significant taxonomy
  #subset_tax <- subset(glom_tax, id %in% sig_neutral_tax$id) 
  
  psuedo_table_name <- paste0(biosample,"_",b,"_table.tsv")
  (write.table(neutral_table, file =psuedo_table_name,sep="\t", row.names = FALSE))
  
  psuedo_map_name <- paste0(biosample,"_",b,"_metadata.txt")
  (write.table(glom_metadata, file =psuedo_map_name,sep="\t", row.names = FALSE))
  
  #psuedo_tax_name <- paste0(biosample,"_",b,"_taxonomy.tsv")
  #(write.table(subset_tax, file =psuedo_tax_name,sep="\t", row.names = FALSE))
}

########### Redo for neutral without filtering for  significant padj ########
print(paste("Writing non-significant taxa table"))
## Subset ASV from each model and signifcance 
sig_neutral_ref = select(filter(non_neutral_tax, model == "neutral" & padj >=0.05),id,p_abundance,freq,padj,Phylum,Class,Order,Family,model)

#subset otu table based on neutral results and rename column names
neutral_table_test <-subset(glom_table2, id  %in%  sig_neutral_ref$id)
la = length(neutral_table_test)
colnames(neutral_table_test)[2:la] <- paste(colnames(neutral_table_test)[2:la], "neutral", sep = "_")
colnames(neutral_table_test)[1] <- "id"

## Create new columns for the metadata data
glom_mapping %>% mutate(psuedo_sample_name = paste(colnames(neutral_table_test)[2:la]),
                            compartment_neutrality = paste(biosample,"neutral", sep = '_'),
                            compartment = paste (biosample),
                            neutrality = paste("neutral")) %>% 
  select(psuedo_sample_name,sample_name_backup, host_scientific_name,expedition_number,political_area,ocean_area,compartment_neutrality,compartment,neutrality) %>% 
  as.data.frame() -> neutral_glom_metadata

## Subset significant results with taxonomy to retain only significant taxonomy
#subset_tax <- subset(glom_tax, id %in% sig_neutral_tax$id) 

## Output the tables for downstream analysis 
psuedo_table_name <- paste0(biosample,"_","neutral_table.tsv")
(write.table(neutral_table_test, file =psuedo_table_name, sep="\t",row.names = FALSE))

psuedo_map_name <- paste0(biosample,"_","neutral_metadata.txt")  ## Remember to combined mapping files in Excel
(write.table(neutral_glom_metadata, file =psuedo_map_name,sep="\t", row.names = FALSE))

################### Merge OTU Tables, Taxonomy & Metadata across neutrality ####################################

############### Import Tables ###############
above_import = paste0(biosample,"_above_table.tsv")
below_import = paste0(biosample,"_below_table.tsv")
neutral_import = paste0(biosample,"_neutral_table.tsv")

above_table <- read.table(above_import,header = TRUE,sep="\t", check.names=FALSE)
neutral_table <- read.table(neutral_import,header = TRUE,sep="\t", check.names=FALSE)
below_table <- read.table(below_import,header = TRUE,sep="\t", check.names=FALSE)

## 1) Merge Tables first 
full_join(above_table, below_table, by ="id") -> merge_table

#merge_table %>% mutate(id = coalesce(id.x,id.y)) %>% 
#  relocate(id) %>%  
#  select(!c(id.x,id.y,Row.names))-> merge_test
merge_table[is.na(merge_table)] <-0

full_join(merge_table,neutral_table, by ="id") -> full_table
#full_table %>% mutate(id = coalesce(id.x,id.y)) %>% 
#  relocate(id) %>%  
#  select(!c(id.x,id.y,Row.names))-> full_table
full_table[is.na(full_table)] <-0

colnames(full_table)[1] <- ""

## Print output files
print(paste("Writing Combined Psudo_table across compartments"))
psudo_table_name <- paste0(biosample,"_combined_psudo_table.tsv")
write.table(full_table, file=psudo_table_name, sep="\t",row.names = FALSE)

############### Merge Taxonomy files ###############
tax2_import = paste0("./taxonomy.tsv")
taxonomy2 <-glom_tax <- read.table(tax2_import, sep = "\t",header = TRUE,check.name=FALSE)
colnames(taxonomy2)[1] <-"id"
neutral_tax <-subset(taxonomy2, id  %in%  full_table$id)
colnames(neutral_tax)[1] <-"Feature ID"
# Print output files
print(paste("Writing Psudo_taxonomy across compartments"))
psudo_tax_name <- paste0(biosample,"_psudo_tax.tsv")
write.table(neutral_tax,psudo_tax_name, sep= "\t", row.names = FALSE)


############### Merge Metadata Files ###############
print(paste("Joining",biosample,"mapping files"))

above_map_import = paste0(biosample,"_above_metadata.txt")
below_map_import = paste0(biosample,"_below_metadata.txt")
neutral_map_import = paste0(biosample,"_neutral_metadata.txt")

above_map <- read.table(above_map_import,header = TRUE,sep="\t", check.names=FALSE)
below_map <- read.table(below_map_import,header = TRUE,sep="\t", check.names=FALSE)
neutral_map <- read.table(neutral_map_import,header = TRUE,sep="\t", check.names=FALSE)

rbind(above_map,below_map,neutral_map) -> combined_map

# Print output files
print(paste("Writing Combined Metadata Filee"))
comb_map_name <- paste0(biosample,"_combined_metadata.txt")
write.table(combined_map,comb_map_name, sep= "\t", row.names = FALSE)


#################### plot neuModel ################################
p_neuM_coral=ggplot(neutral_mod1_w, aes(x=log10(p_abundance), y=freq, color = model))+
  geom_point(shape=19, alpha=0.7, size=2.5)+
  geom_line(aes(x=log10(p_abundance), y=freq.pred),
            color="blue",linetype="solid",linewidth=1)+
  geom_line(aes(x=log10(p_abundance), y=pred.lwr),
            color="blue",linetype="dashed",linewidth=1)+
  geom_line(aes(x=log10(p_abundance), y=pred.upr),
            color="blue",linetype="dashed",linewidth=1)+
  scale_colour_manual(name="OTUs",
                      values=c("#E41A1C", "#999999", "#FF7F00"),
                      breaks=c("above", "neutral", "below"),
                      labels=c("Over-represented", "Neutrally-distributed", "Under-represented"))+
  xlab("Log10 (Mean Relative abundance)") + ylab("Occurance frequency")+
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=12), 
        legend.title = element_text(size=11.5, face = "bold"),
        legend.text = element_text(size=11.5),
        legend.position="none")

r <- round(neutral_mod_w$Rsqr,3)
m_var <- round(neutral_mod_w$m,3)
lb1 <- paste0(r,"\n R2 Value")
lb2 <- paste0(m_var,"\n migration rate (m):")
lb3 <- paste0(biosample,"\n Neutral_Model")

p_neuM_coral1 = p_neuM_coral + annotate("text", x = -3.5, y=0.7, size=5.5, label=lb1, 
                                        fontface="bold", color="black",
                                        parse=FALSE) +
  annotate("text", x = -3.5, y=0.5, size=4.5, label=lb2, 
           fontface="bold", color="black",
           parse=FALSE) +
  annotate("text", x = -3.5, y=0.9, size=4.5, label=lb3, 
           fontface="bold", color="black",
           parse=FALSE)

print(paste("Creating Neutral Figure pdf files"))
## Print Neutral Model Plot
neutral_file_name <- paste0(biosample,"_neutral_model_Plot.pdf")
ggsave(p_neuM_coral1, filename=neutral_file_name)

##################### Generate taxonomic bar graph #########################################

many_col <- c("#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
              "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3","#666666" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8",
              "#4DAF4A" ,"#FF7F00","#984EA3", "#FFFF33", "#A65628","#1B9E77", "#F781BF" ,"#D53E4F","#F46D43", "#FDAE61", "#FEE08B",
              "#E6F598","#ABDDA4","#3288BD","cyan1","#66C2A5","#E69F00", "#56B4E9","#999999","#F0E442","#0072B2","#D55E00","red3",
              "yellow","skyblue","cyan3","deeppink","coral4","#A6CEE3","darkkhaki","brown1","chocolate","darkorchid2","#FF7F00",
              "#66A61E","#FDBF6F", "#FB9A99", "#CAB2D6",  "#B15928","#1B9E77","#D95F02", "#FFFF99","#E7298A", "#7570B3","#666666",
              "#E6AB02", "brown3","#A6761D", "#377EB8", "#4DAF4A" ,"#984EA3","#FFFF33","#FF7F00","black")
tax_rank <- nonneutral_tax %>%
  mutate(model = factor(model, 
                        levels=c("above","neutral","below")),
         Phylum = factor(Phylum),
         Phylum = fct_reorder(Phylum, p_abundance, .desc=TRUE))

### Dotplot Figures ############
print(paste("Creating Non-Neutral Taxa Dotplot" ))
neutral_dot <- ggplot(data = tax_rank, aes(x = model, y = Phylum, size = p_abundance, color = Phylum)) +
  geom_point() + scale_color_manual(values = c(many_col)) + theme_minimal() +
  labs( x="Neutrality", y="Phylum") + guides(color="none")+ ggtitle("Non-Neutral Proteobacteria Bargraph") +
  theme(plot.title = element_text(hjust = 0.5))


print(paste("Creating Neutral Bargraph output files"))
neutral_dotplot_name <- paste0(biosample,"_neutral_dotplot.pdf")
ggsave(neutral_dot, filename=neutral_dotplot_name, width = 20, height = 10, dpi= 600)

print(paste("Finished!"))
