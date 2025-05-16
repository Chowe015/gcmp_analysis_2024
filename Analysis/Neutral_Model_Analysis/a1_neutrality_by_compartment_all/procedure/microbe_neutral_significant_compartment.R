#Loading Required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(minpack.lm)
library(Hmisc)
library(stats4)

sink("PIC_results_log.txt",append=FALSE,split=TRUE)

#Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

glom_table_path <-args[1]
taxonomy_path <- args[2]
biosample <- args[3]

glom_table <- read.csv(glom_table_path, row.names=1, check.names=FALSE)
#print(paste("Imported GCMP OTU Table",otu_table))
glom_taxonomy_table <- read.csv(taxonomy_path, check.name=FALSE)

print(paste("Create Neutral Model Function"))
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  
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
  print(paste("define errors from start"))
  
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

print(paste("Run Neutral Model Function on OTU Table"))

otu_frame = as.data.frame(t(glom_table))

## Neutral fit
neutral_mod_w = sncm.fit(otu_frame)
neutral_mod1_w = sncm.fit(otu_frame,stats = F)

neutral_mod1_w %>% mutate(ci.mean = (pred.upr-pred.lwr)/2,
                        SE = (ci.mean/(1.96)),
                        zstat = (freq/SE),
                        p_value = 2*pnorm(-abs(zstat)),
                        padj = p.adjust(p_value, method = "bonferroni", n=length(neutral_mod1_w$freq))) %>% 
                        as.data.frame() %>%  rownames_to_column("id")-> neutral_results_padj


print(paste("adding model id to asv"))
##Identify above, below and neutral OTUs 
addf = function(freq, pred.upr, pred.lwr, padj) {
  if(freq > pred.upr){
    model = "above"
  }
  else if(freq < pred.lwr){
    model = "below"
  }
  else if(pred.upr >= freq & freq >= pred.lwr & padj < 0.05){
    model = "neutral_significant"
  }
  else if(pred.upr >= freq & freq >= pred.lwr){
    model = "neutral"
  }
  return(model)
}

neutral_results_padj$model = mapply(addf, neutral_results_padj$freq, neutral_results_padj$pred.upr, neutral_results_padj$pred.lwr, neutral_results_padj$padj)  

length(which(neutral_results_padj$model == "above"))
length(which(neutral_results_padj$model == "below"))
length(which(neutral_results_padj$model == "neutral"))
length(which(neutral_results_padj$model == "neutral_significant"))

neutral_results_padj$model <- factor(neutral_results_padj$model, 
                                     levels = c("above", "neutral", "below","neutral_significant"))

print(paste("joining neutral table and ref taxonomy"))
## Inner_join feature table with neutral model families
inner_join(neutral_results_padj,glom_taxonomy_table, by = "id") %>% as.data.frame() -> nonneutral_tax

# Print output files
print(paste("Writing Neutral Model csv"))
neutral_table_name <- paste0(biosample,"_neutral_nonsig_model.csv")
write.csv(neutral_results_padj,neutral_table_name,row.names = FALSE, col.names = TRUE)



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

print(paste("Creating Neutral Results output files"))



## Print Neutral Model Plot
neutral_file_name <- paste0(biosample,"_NeutralModel_nonsig_Plot.pdf")
ggsave(p_neuM_coral1, filename=neutral_file_name, width = 20, height = 10, dpi= 600)
##Generate taxonomic bar graph

many_col <- c("#A6CEE3", "#33A02C","#E7298A","#FDBF6F", "#B2DF8A", "#1F78B4","#E31A1C","#6A3D9A", "#FF7F00", "#FB9A99", "#CAB2D6",
              "#B15928","#1B9E77", "#D95F02", "#FFFF99", "#7570B3","#666666" ,"#66A61E", "#E6AB02","#A6761D","black", "#377EB8",
              "#4DAF4A" ,"#FF7F00","#984EA3", "#FFFF33", "#A65628","#1B9E77", "#F781BF" ,"#D53E4F","#F46D43", "#FDAE61", "#FEE08B",
              "#E6F598","#ABDDA4","#3288BD","cyan1","#66C2A5","#E69F00", "#56B4E9","#999999","#F0E442","#0072B2","#D55E00","red3",
              "yellow","skyblue","cyan3","deeppink","coral4","#A6CEE3","darkkhaki","brown1","chocolate","darkorchid2","#FF7F00",
              "#66A61E","#FDBF6F", "#FB9A99", "#CAB2D6",  "#B15928","#1B9E77","#D95F02", "#FFFF99","#E7298A", "#7570B3","#666666",
              "#E6AB02", "brown3","#A6761D", "#377EB8", "#4DAF4A" ,"#984EA3","#FFFF33","#FF7F00","black")
tax_rank <- nonneutral_tax %>%
  mutate(model = factor(model, 
                        levels=c("above","neutral","neutral_significant","below")),
         Phylum = factor(Phylum),
         Phylum = fct_reorder(Phylum, p_abundance, .desc=TRUE))

neutral_sig_bar <- ggplot(data = tax_rank, aes(x = model, y = p_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position="stack")  + theme_classic() +
  labs( x=" Phylum", y="average relative abundance (p)") +  scale_fill_manual(name=NULL, values = c(many_col))

neutral_bar_name <- paste0(biosample,"_neutral_sig_barplot.pdf")
ggsave(neutral_sig_bar, filename=neutral_bar_name, width = 20, height = 10, dpi= 600)

print(paste("Finished!"))

