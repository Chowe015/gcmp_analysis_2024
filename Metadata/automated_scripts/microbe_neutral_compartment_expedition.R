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

otu_table_path <-args[1]
taxonomy_path <- args[2]
expedition <- args[3]
biosample <- args[4]

otu_table <- read.csv(otu_table_path, row.names=1, check.names=FALSE)
#print(paste("Imported GCMP OTU Table",otu_table))
taxonomy_table <- read.csv(taxonomy_path, row.names=1, check.name=FALSE)
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
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
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

addf = function(freq, pred.upr, pred.lwr) {
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

print(paste("Run Nuetral Model Function on OTU Table"))
## pull in otu table 

otu_frame = as.data.frame(t(otu_table))

## Neutral fit
neutral_mod_w = sncm.fit(otu_frame)
neutral_mod1_w = sncm.fit(otu_frame,stats = F)

#min(sample_sums(glom_col_orb))
print(length(neutral_mod_w$freq))

neutral_mod1_w$model = mapply(addf, neutral_mod1_w$freq, neutral_mod1_w$pred.upr, neutral_mod1_w$pred.lwr)  
length(which(neutral_mod1_w$model == "above"))
length(which(neutral_mod1_w$model == "below"))
length(which(neutral_mod1_w$model == "neutral"))

neutral_mod1_w %>% mutate(ci.mean = (pred.upr-pred.lwr),
                        SE = (ci.mean/(2*1.96)),
                        zstat = (freq/SE),
												p_value = (exp(-0.717*zstat-0.416*zstat^2)),
                        padj = (exp(-0.717*zstat-0.416*zstat^2))*length(neutral_mod1_w$model)) -> neutral_results_padj

print(paste("Print Taxonomy Table from Neutral Model"))

## Assign Taxonomy to table of significant microbes
env_taxa = cbind(as(neutral_results_padj, "data.frame"), as(taxonomy_table[rownames(taxonomy_table), ], "matrix"))
#coral_taxa

#plot neuModel
p_neuM_coral=ggplot(neutral_mod1_w, aes(x=log10(p), y=freq, color = model))+
  geom_point(shape=19, alpha=0.7, size=2.5)+
  geom_line(aes(x=log10(p), y=freq.pred),
            color="blue",linetype="solid",linewidth=1)+
  geom_line(aes(x=log10(p), y=pred.lwr),
            color="blue",linetype="dashed",linewidth=1)+
  geom_line(aes(x=log10(p), y=pred.upr),
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

r <- neutral_mod_w$Rsqr
m_var <- neutral_mod_w$m
lb1 <- paste0(r,"\n R2 Value")
lb2 <- paste0(m_var,"\n migration rate (m):")
lb3 <- paste0(expedition,"_",biosample,"\n Neutral_Model")

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

# Print output files
neutral_table_name <- paste0(expedition,"_",biosample,"_Neutralmodel.csv")
write.csv(env_taxa,neutral_table_name,row.names = TRUE)

## Print Neutral Model Plot
neutral_file_name <- paste0(expedition,"_",biosample,"_NeutralModel_Plot.pdf")
ggsave(p_neuM_coral1, filename=neutral_file_name)

print(paste("Finished!"))

