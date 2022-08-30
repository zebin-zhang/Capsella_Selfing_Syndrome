###############################################################
#        piN/piS       
# The pN/pS ratio is defined as the ratio of the proportion of nonsynonymous polymorphisms to 
# the proportion of synonymous polymorphisms given available sites. 
###############################################################
# install.packages("boot",dep=TRUE)
library(boot)
library(stringr)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
# CG_F
pnps_cg_f <- read.table("InputData/dNdSpiNpiS_CG_F.txt", header = T, sep = "\t", row.names = 1)
nna <- pnps_cg_f$CG_piS != 0
pnps_cg_f <- pnps_cg_f[nna,] # remove genes with piS = 0.
head(pnps_cg_f)

pnps_cg_f_1 <- pnps_cg_f[, c("CG_nb_complete_site", "CG_piN", "CG_piS")]
head(pnps_cg_f_1)
# overlap with all genes in RNAseq
all_F <- read.table("InputData/TMM_Diploids_F.txt", header = T, sep = "\t")
row.names(all_F) <- gsub(".g", "", row.names(all_F))
CGColumns <- str_detect(names(all_F), "CG") # detect which columns contain the "CG" 
CRColumns <- str_detect(names(all_F), "CR") # detect which columns contain the "CR" 
COColumns <- str_detect(names(all_F), "CO") # detect which columns contain the "CO" 

all_F$MeanCG <- rowMeans(all_F[, CGColumns], na.rm = T)
all_F$MeanCR <- rowMeans(all_F[, CRColumns], na.rm = T)
all_F$MeanCO <- rowMeans(all_F[, COColumns], na.rm = T)
head(all_F)
# D_CR & D_CO
Dindex_F <- read.table("OutputData/Flower_D_indcies_all_genes.txt",header = T, sep = "\t")
row.names(Dindex_F) <- gsub(".g", "", row.names(Dindex_F))
head(Dindex_F)

all_F <- merge(all_F, Dindex_F, by=0)
row.names(all_F) <- all_F$Row.names
all_F$Row.names <- NULL
head(all_F)

all_pnps_cg_f <- merge(all_F, pnps_cg_f_1, by = 0) # by row names
row.names(all_pnps_cg_f) <- all_pnps_cg_f$Row.names
all_pnps_cg_f$Row.names <- NULL
head(all_pnps_cg_f)


all_pnps_cg_f <- all_pnps_cg_f[,c("Tissue", "MeanCG", "MeanCR", "MeanCO", "CG_nb_complete_site", "CG_piN", "CG_piS", "D_CR", "D_CO")]
names(all_pnps_cg_f) <-         c("Tissue","MeanCG", "MeanCR", "MeanCO","nb_complete_site", "piN", "piS", "D_CR", "D_CO")
all_pnps_cg_f$Species <- "CG"
all_pnps_cg_f$Tissue <- "Flower"
all_pnps_cg_f$ncs_piN <- all_pnps_cg_f$nb_complete_site * all_pnps_cg_f$piN
all_pnps_cg_f$ncs_piS <- all_pnps_cg_f$nb_complete_site * all_pnps_cg_f$piS

write.table(all_pnps_cg_f, file = "pNpS_data/demo_all_pnps_CG_F.txt", append = F, quote = F, sep = "\t", row.names = T)
# Test whether correlation between expression value and pNpS.
# Slvain saids: As I’ve mentioned, it is also important to control for the level of expression as there usually
# is a negative correlation between expression and piN/piS.
all_pnps_cg_f$ncs_piNpiS <- all_pnps_cg_f$ncs_piN / all_pnps_cg_f$ncs_piS
head(all_pnps_cg_f)
all_pnps_cg_f$series <- NA
all_pnps_cg_f[all_pnps_cg_f$MeanCG < 100, ]$series <- "1-100"
all_pnps_cg_f[all_pnps_cg_f$MeanCG > 100 & all_pnps_cg_f$MeanCG < 500, ]$series <- "100-500"
all_pnps_cg_f[all_pnps_cg_f$MeanCG > 500 & all_pnps_cg_f$MeanCG < 1000, ]$series <- "500-1000"
all_pnps_cg_f[all_pnps_cg_f$MeanCG > 1000 & all_pnps_cg_f$MeanCG < 5000, ]$series <- "1000-5000"
all_pnps_cg_f[all_pnps_cg_f$MeanCG > 5000 , ]$series <- "5000+"

all_pnps_cg_f$series <- factor(all_pnps_cg_f$series, levels = c("1-100", "100-500", "500-1000", "1000-5000", "5000+"))  

ggplot(all_pnps_cg_f, aes(x=MeanCG, y=ncs_piNpiS))+ 
  geom_point(size = 0.5) +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor() +
  labs(x = "TMM expression value", y = "piN / piS") +
  facet_grid(~series, scales="free")
#### seems correlation very low.

# When you consider each gene separately you have a huge variance (because most genes have only a few SNPs).
# You compute the ratio mean(piN) / mean(piS) for the categories you want to compare:
# pN_pS = Sum (nb_complete_site x piN) / Sum (nb_complete_site x piS)

pN_pS_CG_F_all <- sum(all_pnps_cg_f$ncs_piN) / sum(all_pnps_cg_f$ncs_piS )

# CO_F
pnps_co_f <- read.table("InputData/dNdSpiNpiS_CO_F.txt", header = T, sep = "\t", row.names = 1)
nna <- pnps_co_f$CO_piS != 0
pnps_co_f <- pnps_co_f[nna,] # remove genes with piS = 0.
head(pnps_co_f)

pnps_co_f_1 <- pnps_co_f[, c("CO_nb_complete_site", "CO_piN", "CO_piS")]
head(pnps_co_f_1)
# overlap with all genes in RNAseq
head(all_F)
all_pnps_co_f <- merge(all_F, pnps_co_f_1, by = 0) # by row names
row.names(all_pnps_co_f) <- all_pnps_co_f$Row.names
all_pnps_co_f$Row.names <- NULL

all_pnps_co_f <- all_pnps_co_f[,c("Tissue", "MeanCG", "MeanCR", "MeanCO","CO_nb_complete_site", "CO_piN", "CO_piS", "D_CR", "D_CO")]
names(all_pnps_co_f) <- c("Tissue", "MeanCG", "MeanCR", "MeanCO","nb_complete_site", "piN", "piS", "D_CR", "D_CO")
all_pnps_co_f$Species <- "CO"
all_pnps_co_f$Tissue <- "Flower"
all_pnps_co_f$ncs_piN <- all_pnps_co_f$nb_complete_site * all_pnps_co_f$piN
all_pnps_co_f$ncs_piS <- all_pnps_co_f$nb_complete_site * all_pnps_co_f$piS
head(all_pnps_co_f)
# When you consider each gene separately you have a huge variance (because most genes have only a few SNPs).
# You compute the ratio mean(piN) / mean(piS) for the categories you want to compare:
# pN_pS = Sum (nb_complete_site x piN) / Sum (nb_complete_site x piS)

pN_pS_CO_F_all <- sum(all_pnps_co_f$ncs_piN) / sum(all_pnps_co_f$ncs_piS )

# CR_F
pnps_cr_f <- read.table("InputData/dNdSpiNpiS_CR_F.txt", header = T, sep = "\t", row.names = 1)
nna <- pnps_cr_f$CR_piS != 0
pnps_cr_f <- pnps_cr_f[nna,] # remove genes with piS = 0.
head(pnps_cr_f)

pnps_cr_f_1 <- pnps_cr_f[, c("CR_nb_complete_site", "CR_piN", "CR_piS")]
head(pnps_cr_f_1)
# overlap with all genes in RNAseq
all_pnps_cr_f <- merge(all_F, pnps_cr_f_1, by = 0) # by row names
row.names(all_pnps_cr_f) <- all_pnps_cr_f$Row.names
all_pnps_cr_f$Row.names <- NULL

all_pnps_cr_f <- all_pnps_cr_f[,c(c("Tissue", "MeanCG", "MeanCR", "MeanCO","CR_nb_complete_site", "CR_piN", "CR_piS", "D_CR", "D_CO"))]
names(all_pnps_cr_f) <- c(c("Tissue", "MeanCG", "MeanCR", "MeanCO","nb_complete_site", "piN", "piS", "D_CR", "D_CO"))
all_pnps_cr_f$Species <- "CR"
all_pnps_cr_f$Tissue <- "Flower"
all_pnps_cr_f$ncs_piN <- all_pnps_cr_f$nb_complete_site * all_pnps_cr_f$piN
all_pnps_cr_f$ncs_piS <- all_pnps_cr_f$nb_complete_site * all_pnps_cr_f$piS
head(all_pnps_cr_f)
# When you consider each gene separately you have a huge variance (because most genes have only a few SNPs).
# You compute the ratio mean(piN) / mean(piS) for the categories you want to compare:
# pN_pS = Sum (nb_complete_site x piN) / Sum (nb_complete_site x piS)

pN_pS_CR_F_all <- sum(all_pnps_cr_f$ncs_piN) / sum(all_pnps_cr_f$ncs_piS )

# bootstrap CI 95%

pnps <- function(d, i){
  d2 <- d[i,]
  pn <- d2$nb_complete_site * d2$piN
  ps <- d2$nb_complete_site * d2$piS
  return(sum(pn) / sum(ps))
}

bootstrap_all_pnps_cg_f <- boot(all_pnps_cg_f, pnps, R=1000)
summary(bootstrap_all_pnps_cg_f)
boot.ci(boot.out=bootstrap_all_pnps_cg_f,type=c("norm","basic","perc","bca"))
sort(bootstrap_all_pnps_cg_f$t)[c(25,975)]

bootstrap_all_pnps_co_f <- boot(all_pnps_co_f, pnps, R=10000)
summary(bootstrap_all_pnps_co_f)
sort(bootstrap_all_pnps_co_f$t)[c(250,9750)]
boot.ci(boot.out=bootstrap_all_pnps_co_f,type=c("norm","basic","perc","bca"))

bootstrap_all_pnps_cr_f <- boot(all_pnps_cr_f, pnps, R=10000)
summary(bootstrap_all_pnps_cr_f)
sort(bootstrap_all_pnps_cr_f$t)[c(250,9750)]
boot.ci(boot.out=bootstrap_all_pnps_cr_f,type=c("norm","basic","perc","bca"))

## piN/piS and coDEGs
## CG Overlap with coDEGs_CRCO
head(all_pnps_cg_f)
dim(all_pnps_cg_f)
coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_F) <- gsub(".g", "", row.names(coDEG_F))
coDEG_F$coDEGs <- "coDEGs"
coDEG_F$Tissue <- NULL
head(coDEG_F)
codeg_pnps_cg_f <- merge(all_pnps_cg_f, coDEG_F, by = 0, all = T)
codeg_pnps_cg_f[is.na(codeg_pnps_cg_f$DE_1),]$coDEGs <- "Others"
codeg_pnps_cg_f[is.na(codeg_pnps_cg_f$DE_1),]$DE_1 <- "Others"
codeg_pnps_cg_f <- codeg_pnps_cg_f[!is.na(codeg_pnps_cg_f$piS),]
codeg_pnps_cg_f$Tissue <- "Flower"
head(codeg_pnps_cg_f)
dim(codeg_pnps_cg_f)
row_codeg_f <- codeg_pnps_cg_f$coDEGs == "coDEGs"
codeg_down_cg_f <- codeg_pnps_cg_f$DE_1 == "Down"
codeg_up_cg_f <- codeg_pnps_cg_f$DE_1 == "Up"
row_other_f <- codeg_pnps_cg_f$coDEGs == "Others"
d_pos_cg_f <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0 & codeg_pnps_cg_f$D_CO > 0 # D_co D_cr positive overlap
d_neg_cg_f <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < 0 & codeg_pnps_cg_f$D_CO < 0 # D_co D_cr negative overlap
d_neg_cg_f_0.1 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.1 & codeg_pnps_cg_f$D_CO < -0.1 # D_co D_cr negative overlap
d_neg_cg_f_0.3 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.3 & codeg_pnps_cg_f$D_CO < -0.3 # D_co D_cr negative overlap
d_neg_cg_f_0.5 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.5 & codeg_pnps_cg_f$D_CO < -0.5 # D_co D_cr negative overlap
d_neg_cg_f_0.7 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.7 & codeg_pnps_cg_f$D_CO < -0.7 # D_co D_cr negative overlap
d_neg_cg_f_0.8 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.8 & codeg_pnps_cg_f$D_CO < -0.8 # D_co D_cr negative overlap
d_neg_cg_f_0.9 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR < -0.9 & codeg_pnps_cg_f$D_CO < -0.9 # D_co D_cr negative overlap
d_pos_cg_f_0.1 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.1 & codeg_pnps_cg_f$D_CO > 0.1 # D_co D_cr negative overlap
d_pos_cg_f_0.3 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.3 & codeg_pnps_cg_f$D_CO > 0.3 # D_co D_cr negative overlap
d_pos_cg_f_0.5 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.5 & codeg_pnps_cg_f$D_CO > 0.5 # D_co D_cr negative overlap
d_pos_cg_f_0.7 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.7 & codeg_pnps_cg_f$D_CO > 0.7 # D_co D_cr negative overlap
d_pos_cg_f_0.8 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.8 & codeg_pnps_cg_f$D_CO > 0.8 # D_co D_cr negative overlap
d_pos_cg_f_0.9 <- codeg_pnps_cg_f$coDEGs == "coDEGs" & codeg_pnps_cg_f$D_CR > 0.9 & codeg_pnps_cg_f$D_CO > 0.9 # D_co D_cr negative overlap


pN_pS_CG_F_codeg <- sum(codeg_pnps_cg_f[row_codeg_f,]$ncs_piN) / sum(codeg_pnps_cg_f[row_codeg_f,]$ncs_piS )
pN_pS_CG_F_codeg_up <- sum(codeg_pnps_cg_f[codeg_up_cg_f,]$ncs_piN) / sum(codeg_pnps_cg_f[codeg_up_cg_f,]$ncs_piS )
pN_pS_CG_F_codeg_down <- sum(codeg_pnps_cg_f[codeg_down_cg_f,]$ncs_piN) / sum(codeg_pnps_cg_f[codeg_down_cg_f,]$ncs_piS )
pN_pS_CG_F_other <- sum(codeg_pnps_cg_f[row_other_f,]$ncs_piN) / sum(codeg_pnps_cg_f[row_other_f,]$ncs_piS )
pN_pS_CG_F_d_pos <- sum(codeg_pnps_cg_f[d_pos_cg_f,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f,]$ncs_piS )
pN_pS_CG_F_d_neg <- sum(codeg_pnps_cg_f[d_neg_cg_f,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f,]$ncs_piS )
pN_pS_CG_F_d_neg_0.1 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.1,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.1,]$ncs_piS )
pN_pS_CG_F_d_neg_0.3 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.3,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.3,]$ncs_piS )
pN_pS_CG_F_d_neg_0.5 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.5,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.5,]$ncs_piS )
pN_pS_CG_F_d_neg_0.7 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.7,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.7,]$ncs_piS )
pN_pS_CG_F_d_neg_0.8 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.8,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.8,]$ncs_piS )
pN_pS_CG_F_d_neg_0.9 <- sum(codeg_pnps_cg_f[d_neg_cg_f_0.9,]$ncs_piN) / sum(codeg_pnps_cg_f[d_neg_cg_f_0.9,]$ncs_piS )
pN_pS_CG_F_d_pos_0.1 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.1,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.1,]$ncs_piS )
pN_pS_CG_F_d_pos_0.3 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.3,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.3,]$ncs_piS )
pN_pS_CG_F_d_pos_0.5 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.5,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.5,]$ncs_piS )
pN_pS_CG_F_d_pos_0.7 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.7,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.7,]$ncs_piS )
pN_pS_CG_F_d_pos_0.8 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.8,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.8,]$ncs_piS )
pN_pS_CG_F_d_pos_0.9 <- sum(codeg_pnps_cg_f[d_pos_cg_f_0.9,]$ncs_piN) / sum(codeg_pnps_cg_f[d_pos_cg_f_0.9,]$ncs_piS )

## CR Overlap with coDEGs_CRCO
head(all_pnps_cr_f)
# coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
# row.names(coDEG_F) <- gsub(".g", "", row.names(coDEG_F))
# coDEG_F$coDEGs <- "coDEGs"
# coDEG_F$Tissue <- NULL
head(coDEG_F)
codeg_pnps_cr_f <- merge(all_pnps_cr_f, coDEG_F, by = 0, all = T)
codeg_pnps_cr_f[is.na(codeg_pnps_cr_f$DE_1),]$coDEGs <- "Others"
codeg_pnps_cr_f[is.na(codeg_pnps_cr_f$DE_1),]$DE_1 <- "Others"
codeg_pnps_cr_f <- codeg_pnps_cr_f[!is.na(codeg_pnps_cr_f$piS),]
codeg_pnps_cr_f$Tissue <- "Flower"
codeg_pnps_cr_f_backup <- codeg_pnps_cr_f
codeg_pnps_cr_f <- codeg_pnps_cr_f[codeg_pnps_cr_f$MeanCR < 5000,]
head(codeg_pnps_cr_f)
row_codeg_cr_f <- codeg_pnps_cr_f$coDEGs == "coDEGs"
codeg_up_cr_f <- codeg_pnps_cr_f$DE_1 == "Up"
codeg_down_cr_f <- codeg_pnps_cr_f$DE_1 == "Down"
row_other_cr_f <- codeg_pnps_cr_f$coDEGs == "Others"
d_pos_cr_f <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0 & codeg_pnps_cr_f$D_CO > 0 # D_co D_cr positive overlap
d_neg_cr_f <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < 0 & codeg_pnps_cr_f$D_CO < 0 # D_co D_cr negative overlap
d_neg_cr_f_0.1 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.1 & codeg_pnps_cr_f$D_CO < -0.1 # D_co D_cr negative overlap
d_neg_cr_f_0.3 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.3 & codeg_pnps_cr_f$D_CO < -0.3 # D_co D_cr negative overlap
d_neg_cr_f_0.5 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.5 & codeg_pnps_cr_f$D_CO < -0.5 # D_co D_cr negative overlap
d_neg_cr_f_0.7 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.7 & codeg_pnps_cr_f$D_CO < -0.7 # D_co D_cr negative overlap
d_neg_cr_f_0.8 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.8 & codeg_pnps_cr_f$D_CO < -0.8 # D_co D_cr negative overlap
d_neg_cr_f_0.9 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR < -0.9 & codeg_pnps_cr_f$D_CO < -0.9 # D_co D_cr negative overlap
d_pos_cr_f_0.1 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.1 & codeg_pnps_cr_f$D_CO  > 0.1 # D_co D_cr negative overlap
d_pos_cr_f_0.3 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.3 & codeg_pnps_cr_f$D_CO  > 0.3 # D_co D_cr negative overlap
d_pos_cr_f_0.5 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.5 & codeg_pnps_cr_f$D_CO  > 0.5 # D_co D_cr negative overlap
d_pos_cr_f_0.7 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.7 & codeg_pnps_cr_f$D_CO  > 0.7 # D_co D_cr negative overlap
d_pos_cr_f_0.8 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.8 & codeg_pnps_cr_f$D_CO  > 0.8 # D_co D_cr negative overlap
d_pos_cr_f_0.9 <- codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$D_CR > 0.9 & codeg_pnps_cr_f$D_CO  > 0.9 # D_co D_cr negative overlap

pN_pS_CR_F_codeg <- sum(codeg_pnps_cr_f[row_codeg_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[row_codeg_cr_f,]$ncs_piS )
pN_pS_CR_F_codeg_up <- sum(codeg_pnps_cr_f[codeg_up_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[codeg_up_cr_f,]$ncs_piS )
pN_pS_CR_F_codeg_down <- sum(codeg_pnps_cr_f[codeg_down_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[codeg_down_cr_f,]$ncs_piS )
pN_pS_CR_F_other <- sum(codeg_pnps_cr_f[row_other_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[row_other_cr_f,]$ncs_piS )
pN_pS_CR_F_d_pos <- sum(codeg_pnps_cr_f[d_pos_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f,]$ncs_piS )
pN_pS_CR_F_d_neg <- sum(codeg_pnps_cr_f[d_neg_cr_f,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f,]$ncs_piS )
pN_pS_CR_F_d_neg_0.1 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.1,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.1,]$ncs_piS )
pN_pS_CR_F_d_neg_0.3 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.3,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.3,]$ncs_piS )
pN_pS_CR_F_d_neg_0.5 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.5,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.5,]$ncs_piS )
pN_pS_CR_F_d_neg_0.7 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.7,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.7,]$ncs_piS )
pN_pS_CR_F_d_neg_0.8 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.8,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.8,]$ncs_piS )
pN_pS_CR_F_d_neg_0.9 <- sum(codeg_pnps_cr_f[d_neg_cr_f_0.9,]$ncs_piN) / sum(codeg_pnps_cr_f[d_neg_cr_f_0.9,]$ncs_piS )
pN_pS_CR_F_d_pos_0.1 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.1,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.1,]$ncs_piS )
pN_pS_CR_F_d_pos_0.3 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.3,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.3,]$ncs_piS )
pN_pS_CR_F_d_pos_0.5 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.5,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.5,]$ncs_piS )
pN_pS_CR_F_d_pos_0.7 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.7,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.7,]$ncs_piS )
pN_pS_CR_F_d_pos_0.8 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.8,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.8,]$ncs_piS )
pN_pS_CR_F_d_pos_0.9 <- sum(codeg_pnps_cr_f[d_pos_cr_f_0.9,]$ncs_piN) / sum(codeg_pnps_cr_f[d_pos_cr_f_0.9,]$ncs_piS )

testdd <- codeg_pnps_cr_f[d_neg_cr_f_0.9,]
## CO Overlap with coDEGs_CRCO
head(all_pnps_co_f)
# coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
# row.names(coDEG_F) <- gsub(".g", "", row.names(coDEG_F))
# coDEG_F$coDEGs <- "coDEGs"
# coDEG_F$Tissue <- NULL
head(coDEG_F)
codeg_pnps_co_f <- merge(all_pnps_co_f, coDEG_F, by = 0, all = T)
codeg_pnps_co_f[is.na(codeg_pnps_co_f$DE_1),]$coDEGs <- "Others"
codeg_pnps_co_f[is.na(codeg_pnps_co_f$DE_1),]$DE_1 <- "Others"
codeg_pnps_co_f <- codeg_pnps_co_f[!is.na(codeg_pnps_co_f$piS),]
codeg_pnps_co_f$Tissue <- "Flower"
head(codeg_pnps_co_f)
row_codeg_co_f <- codeg_pnps_co_f$coDEGs == "coDEGs"
codeg_up_co_f <- codeg_pnps_co_f$DE_1 == "Up"
codeg_down_co_f <- codeg_pnps_co_f$DE_1 == "Down"
row_other_co_f <- codeg_pnps_co_f$coDEGs == "Others"
d_pos_co_f <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0 & codeg_pnps_co_f$D_CO > 0 # D_co D_cr positive overlap
d_neg_co_f <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < 0 & codeg_pnps_co_f$D_CO < 0 # D_co D_cr negative overlap
d_neg_co_f_0.1 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.1 & codeg_pnps_co_f$D_CO < -0.1 # D_co D_cr negative overlap
d_neg_co_f_0.3 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.3 & codeg_pnps_co_f$D_CO < -0.3 # D_co D_cr negative overlap
d_neg_co_f_0.5 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.5 & codeg_pnps_co_f$D_CO < -0.5 # D_co D_cr negative overlap
d_neg_co_f_0.7 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.7 & codeg_pnps_co_f$D_CO < -0.7 # D_co D_cr negative overlap
d_neg_co_f_0.8 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.8 & codeg_pnps_co_f$D_CO < -0.8 # D_co D_cr negative overlap
d_neg_co_f_0.9 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR < -0.9 & codeg_pnps_co_f$D_CO < -0.9 # D_co D_cr negative overlap
d_pos_co_f_0.1 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.1 & codeg_pnps_co_f$D_CO > 0.1 # D_co D_cr negative overlap
d_pos_co_f_0.3 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.3 & codeg_pnps_co_f$D_CO > 0.3 # D_co D_cr negative overlap
d_pos_co_f_0.5 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.5 & codeg_pnps_co_f$D_CO > 0.5 # D_co D_cr negative overlap
d_pos_co_f_0.7 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.7 & codeg_pnps_co_f$D_CO > 0.7 # D_co D_cr negative overlap
d_pos_co_f_0.8 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.8 & codeg_pnps_co_f$D_CO > 0.8 # D_co D_cr negative overlap
d_pos_co_f_0.9 <- codeg_pnps_co_f$coDEGs == "coDEGs" & codeg_pnps_co_f$D_CR > 0.9 & codeg_pnps_co_f$D_CO > 0.9 # D_co D_cr negative overlap


pN_pS_CO_F_codeg <- sum(codeg_pnps_co_f[row_codeg_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[row_codeg_co_f,]$ncs_piS )
pN_pS_CO_F_codeg_up <- sum(codeg_pnps_co_f[codeg_up_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[codeg_up_co_f,]$ncs_piS )
pN_pS_CO_F_codeg_down <- sum(codeg_pnps_co_f[codeg_down_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[codeg_down_co_f,]$ncs_piS )
pN_pS_CO_F_other <- sum(codeg_pnps_co_f[row_other_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[row_other_co_f,]$ncs_piS )
pN_pS_CO_F_d_pos <- sum(codeg_pnps_co_f[d_pos_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f,]$ncs_piS )
pN_pS_CO_F_d_neg <- sum(codeg_pnps_co_f[d_neg_co_f,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f,]$ncs_piS )
pN_pS_CO_F_d_neg_0.1 <- sum(codeg_pnps_co_f[d_neg_co_f_0.1,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.1,]$ncs_piS )
pN_pS_CO_F_d_neg_0.3 <- sum(codeg_pnps_co_f[d_neg_co_f_0.3,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.3,]$ncs_piS )
pN_pS_CO_F_d_neg_0.5 <- sum(codeg_pnps_co_f[d_neg_co_f_0.5,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.5,]$ncs_piS )
pN_pS_CO_F_d_neg_0.7 <- sum(codeg_pnps_co_f[d_neg_co_f_0.7,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.7,]$ncs_piS )
pN_pS_CO_F_d_neg_0.8 <- sum(codeg_pnps_co_f[d_neg_co_f_0.8,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.8,]$ncs_piS )
pN_pS_CO_F_d_neg_0.9 <- sum(codeg_pnps_co_f[d_neg_co_f_0.9,]$ncs_piN) / sum(codeg_pnps_co_f[d_neg_co_f_0.9,]$ncs_piS )
pN_pS_CO_F_d_pos_0.1 <- sum(codeg_pnps_co_f[d_pos_co_f_0.1,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.1,]$ncs_piS )
pN_pS_CO_F_d_pos_0.3 <- sum(codeg_pnps_co_f[d_pos_co_f_0.3,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.3,]$ncs_piS )
pN_pS_CO_F_d_pos_0.5 <- sum(codeg_pnps_co_f[d_pos_co_f_0.5,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.5,]$ncs_piS )
pN_pS_CO_F_d_pos_0.7 <- sum(codeg_pnps_co_f[d_pos_co_f_0.7,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.7,]$ncs_piS )
pN_pS_CO_F_d_pos_0.8 <- sum(codeg_pnps_co_f[d_pos_co_f_0.8,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.8,]$ncs_piS )
pN_pS_CO_F_d_pos_0.9 <- sum(codeg_pnps_co_f[d_pos_co_f_0.9,]$ncs_piN) / sum(codeg_pnps_co_f[d_pos_co_f_0.9,]$ncs_piS )

# CG MST genes
head(all_pnps_cg_f)
MST <- read.table("OutputData/482_MST_genes_overlap_Slotte_Wozniak.txt", header = T, sep = "\t")
MST$gene_tpye <- "MST"
MST <- MST[, c("DE", "gene_tpye")]
row.names(MST) <- gsub(".g", "", row.names(MST))
head(MST)
mst_pnps_cg_f <- merge(all_pnps_cg_f, MST, by = 0, all = T)
mst_pnps_cg_f[is.na(mst_pnps_cg_f$gene_tpye),]$gene_tpye <- "Others"
mst_pnps_cg_f[is.na(mst_pnps_cg_f$DE),]$DE <- "Others"
mst_pnps_cg_f <- mst_pnps_cg_f[!is.na(mst_pnps_cg_f$piS),]
mst_pnps_cg_f$Tissue <- "Flower"
head(mst_pnps_cg_f)
row_mst_cg_f <- mst_pnps_cg_f$gene_tpye == "MST"
row_down_cg_f <- mst_pnps_cg_f$DE == "Up" # convert CG/CR > CR/CG
row_up_cg_f <- mst_pnps_cg_f$DE == "Down" # convert CG/CR
row_oth_cg_f <- mst_pnps_cg_f$gene_tpye == "Others"

pN_pS_CG_F_mst <- sum(mst_pnps_cg_f[row_mst_cg_f,]$ncs_piN) / sum(mst_pnps_cg_f[row_mst_cg_f,]$ncs_piS )
pN_pS_CG_F_up <- sum(mst_pnps_cg_f[row_up_cg_f,]$ncs_piN) / sum(mst_pnps_cg_f[row_up_cg_f,]$ncs_piS)
pN_pS_CG_F_down <- sum(mst_pnps_cg_f[row_down_cg_f,]$ncs_piN) / sum(mst_pnps_cg_f[row_down_cg_f,]$ncs_piS )
pN_pS_CG_F_oth <- sum(mst_pnps_cg_f[row_oth_cg_f,]$ncs_piN) / sum(mst_pnps_cg_f[row_oth_cg_f,]$ncs_piS )

# CR MST genes
head(all_pnps_cr_f)
# MST <- read.table("OutputData/482_MST_genes_overlap_Slotte_Wozniak.txt", header = T, sep = "\t")
# MST$gene_tpye <- "MST"
# MST <- MST[, c("DE", "gene_tpye")]
# row.names(MST) <- gsub(".g", "", row.names(MST))
head(MST)
mst_pnps_cr_f <- merge(all_pnps_cr_f, MST, by = 0, all = T)
mst_pnps_cr_f[is.na(mst_pnps_cr_f$gene_tpye),]$gene_tpye <- "Others"
mst_pnps_cr_f[is.na(mst_pnps_cr_f$DE),]$DE <- "Others"
mst_pnps_cr_f <- mst_pnps_cr_f[!is.na(mst_pnps_cr_f$piS),]
mst_pnps_cr_f$Tissue <- "Flower"
head(mst_pnps_cr_f)
row_mst_cr_f <- mst_pnps_cr_f$gene_tpye == "MST"
row_down_cr_f <- mst_pnps_cr_f$DE == "Up" # convert CG/CR > CR/CG
row_up_cr_f <- mst_pnps_cr_f$DE == "Down" # convert CG/CR > CR/CG
row_oth_cr_f <- mst_pnps_cr_f$gene_tpye == "Others"

pN_pS_CR_F_mst <- sum(mst_pnps_cr_f[row_mst_cr_f,]$ncs_piN) / sum(mst_pnps_cr_f[row_mst_cr_f,]$ncs_piS )
pN_pS_CR_F_up <- sum(mst_pnps_cr_f[row_up_cr_f,]$ncs_piN) / sum(mst_pnps_cr_f[row_up_cr_f,]$ncs_piS)
pN_pS_CR_F_down <- sum(mst_pnps_cr_f[row_down_cr_f,]$ncs_piN) / sum(mst_pnps_cr_f[row_down_cr_f,]$ncs_piS )
pN_pS_CR_F_oth <- sum(mst_pnps_cr_f[row_oth_cr_f,]$ncs_piN) / sum(mst_pnps_cr_f[row_oth_cr_f,]$ncs_piS )

# CO MST genes
head(all_pnps_co_f)
# MST <- read.table("OutputData/482_MST_genes_overlap_Slotte_Wozniak.txt", header = T, sep = "\t")
# MST$gene_tpye <- "MST"
# MST <- MST[, c("DE", "gene_tpye")]
# row.names(MST) <- gsub(".g", "", row.names(MST))
head(MST)
mst_pnps_co_f <- merge(all_pnps_co_f, MST, by = 0, all = T)
mst_pnps_co_f[is.na(mst_pnps_co_f$gene_tpye),]$gene_tpye <- "Others"
mst_pnps_co_f[is.na(mst_pnps_co_f$DE),]$DE <- "Others"
mst_pnps_co_f <- mst_pnps_co_f[!is.na(mst_pnps_co_f$piS),]
mst_pnps_co_f$Tissue <- "Flower"
head(mst_pnps_co_f)
row_mst_co_f <- mst_pnps_co_f$gene_tpye == "MST"
row_down_co_f <- mst_pnps_co_f$DE == "Up" # convert CG/CR > CR/CG
row_up_co_f <- mst_pnps_co_f$DE == "Down" # convert CG/CR > CR/CG
row_oth_co_f <- mst_pnps_co_f$gene_tpye == "Others"

pN_pS_CO_F_mst <- sum(mst_pnps_co_f[row_mst_co_f,]$ncs_piN) / sum(mst_pnps_co_f[row_mst_co_f,]$ncs_piS )
pN_pS_CO_F_up <- sum(mst_pnps_co_f[row_up_co_f,]$ncs_piN) / sum(mst_pnps_co_f[row_up_co_f,]$ncs_piS)
pN_pS_CO_F_down <- sum(mst_pnps_co_f[row_down_co_f,]$ncs_piN) / sum(mst_pnps_co_f[row_down_co_f,]$ncs_piS )
pN_pS_CO_F_oth <- sum(mst_pnps_co_f[row_oth_co_f,]$ncs_piN) / sum(mst_pnps_co_f[row_oth_co_f,]$ncs_piS )





























pnps_cg_f[nna,]$CG_F_pNpS <- pnps_cg_f[nna,]$CG_piN / pnps_cg_f[nna,]$CG_piS

row.names(pnps_cg_f_1) <- pnps_cg_f_1$Contig_name
head(pnps_cg_f_1)
# overlap with all genes in RNAseq
all_F <- read.table("InputData/TMM_Diploids_F.txt", header = T, sep = "\t")
row.names(all_F) <- gsub(".g", "", row.names(all_F))
head(all_F)
all_pnps_cg_f <- merge(all_F, pnps_cg_f_1, by = 0) # by row names
head(all_pnps_cg_f)
row.names(all_pnps_cg_f) <- all_pnps_cg_f$Row.names
all_pnps_cg_f$Row.names <- NULL
dim(all_pnps_cg_f)
all_pnps_cg_f <- all_pnps_cg_f[, c("Contig_name", "CG_piN", "CG_piS", "CG_F_pNpS")]
head(all_pnps_cg_f)


# CO_F
pnps_co_f <- read.table("InputData/dNdSpiNpiS_CO_F.txt", header = T, sep = "\t")
head(pnps_co_f)
pnps_co_f$CO_F_pNpS <- 0
nna <- pnps_co_f$CO_piS != 0
pnps_co_f[nna,]$CO_F_pNpS <- pnps_co_f[nna,]$CO_piN / pnps_co_f[nna,]$CO_piS
pnps_co_f_1 <- pnps_co_f[, c("Contig_name", "CO_piN", "CO_piS", "CO_F_pNpS")]
row.names(pnps_co_f_1) <- pnps_co_f_1$Contig_name
head(pnps_co_f_1)
# overlap with all genes in RNAseq
all_pnps_co_f <- merge(all_F, pnps_co_f_1, by = 0) # by row names
head(all_pnps_co_f)
row.names(all_pnps_co_f) <- all_pnps_co_f$Row.names
all_pnps_co_f$Row.names <- NULL
dim(all_pnps_co_f)
all_pnps_co_f <- all_pnps_co_f[, c("Contig_name", "CO_piN", "CO_piS", "CO_F_pNpS")]
head(all_pnps_co_f)
# Overlap with coDEGs_CRCO
coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_F) <- gsub(".g", "", row.names(coDEG_F))
coDEG_F$coDEGs <- "coDEGs"
head(coDEG_F)
codeg_pnps_co_f <- merge(all_pnps_co_f, coDEG_F, by = 0, all = T)
codeg_pnps_co_f[is.na(codeg_pnps_co_f$DE_1),]$coDEGs <- "Others"
codeg_pnps_co_f$Tissue <- "Flower"
head(codeg_pnps_co_f)
mean(codeg_pnps_co_f[codeg_pnps_co_f$coDEGs == "coDEGs",]$CO_F_pNpS, na.rm = T)
mean(codeg_pnps_co_f[codeg_pnps_co_f$coDEGs == "Others",]$CO_F_pNpS, na.rm = T)

### CR
pnps_cr_f <- read.table("InputData/dNdSpiNpiS_CR_F.txt", header = T, sep = "\t")
head(pnps_cr_f)
pnps_cr_f$CR_F_pNpS <- 0
nna <- pnps_cr_f$CR_piS != 0
pnps_cr_f[nna,]$CR_F_pNpS <- pnps_cr_f[nna,]$CR_piN / pnps_cr_f[nna,]$CR_piS
pnps_cr_f_1 <- pnps_cr_f[, c("Contig_name", "CR_piN", "CR_piS", "CR_F_pNpS")]
row.names(pnps_cr_f_1) <- pnps_cr_f_1$Contig_name
head(pnps_cr_f_1)
# overlap with all genes in RNAseq
all_pnps_cr_f <- merge(all_F, pnps_cr_f_1, by = 0) # by row names
head(all_pnps_cr_f)
row.names(all_pnps_cr_f) <- all_pnps_cr_f$Row.names
all_pnps_cr_f$Row.names <- NULL
dim(all_pnps_cr_f)
all_pnps_cr_f <- all_pnps_cr_f[, c("Contig_name", "CR_piN", "CR_piS", "CR_F_pNpS")]
head(all_pnps_cr_f)
# Overlap with coDEGs_CRCO
coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_F) <- gsub(".g", "", row.names(coDEG_F))
coDEG_F$coDEGs <- "coDEGs"
head(coDEG_F)
codeg_pnps_cr_f <- merge(all_pnps_cr_f, coDEG_F, by = 0, all = T)
codeg_pnps_cr_f[is.na(codeg_pnps_cr_f$DE_1),]$coDEGs <- "Others"
codeg_pnps_cr_f$Tissue <- "Flower"
head(codeg_pnps_cr_f)
median(codeg_pnps_cr_f[codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$CR_F_pNpS != 0,]$CR_F_pNpS, na.rm = T)
median(codeg_pnps_cr_f[codeg_pnps_cr_f$coDEGs == "Others" & codeg_pnps_cr_f$CR_F_pNpS != 0,]$CR_F_pNpS, na.rm = T)
mean(codeg_pnps_cr_f[codeg_pnps_cr_f$coDEGs == "coDEGs" & codeg_pnps_cr_f$CR_F_pNpS != 0,]$CR_F_pNpS, na.rm = T)
mean(codeg_pnps_cr_f[codeg_pnps_cr_f$coDEGs == "Others" & codeg_pnps_cr_f$CR_F_pNpS != 0,]$CR_F_pNpS, na.rm = T)

# combine CG, CO, CR
head(codeg_pnps_cg_f)
cg_pnps_f <- codeg_pnps_cg_f[,c("Contig_name","CG_piN", "CG_piS", "CG_F_pNpS", "Tissue", "coDEGs")]
cg_pnps_f$Species <- "CG"
cg_pnps_f <- cg_pnps_f[!(is.na(cg_pnps_f$CG_F_pNpS)),]
names(cg_pnps_f) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(cg_pnps_f)

co_pnps_f <- codeg_pnps_co_f[,c("Contig_name","CO_piN", "CO_piS", "CO_F_pNpS", "Tissue", "coDEGs")]
co_pnps_f$Species <- "CO"
co_pnps_f <- co_pnps_f[!(is.na(co_pnps_f$CO_F_pNpS)),]
names(co_pnps_f) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(co_pnps_f)

cr_pnps_f <- codeg_pnps_cr_f[,c("Contig_name","CR_piN", "CR_piS", "CR_F_pNpS", "Tissue", "coDEGs")]
cr_pnps_f$Species <- "CR"
cr_pnps_f <- cr_pnps_f[!(is.na(cr_pnps_f$CR_F_pNpS)),]
names(cr_pnps_f) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(cr_pnps_f)

pnps_f <- rbind(cg_pnps_f, co_pnps_f, cr_pnps_f)
head(pnps_f)

my_comparisons <- list( c("coDEGs", "Others"))

ggplot(pnps_f[pnps_f$piS != 0,], aes(x=coDEGs, y=piNpiS, color=coDEGs)) + 
  geom_boxplot() +
  scale_color_manual(values=c("blue","orange")) +
  facet_grid( ~ Species) +
  stat_compare_means(label.y = 2.5) +
  stat_summary(fun = mean, geom = "point", shape=5, size=1.5, color = "black") +
  #stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
  ylim(0,2.5) +
  theme_bw() +
  labs(x = "Flowers", y = "piN / piS")

ggplot(pnps_f[pnps_f$piNpiS != 0,], aes(x=piNpiS, fill=coDEGs)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +
  xlim(0,2.5) +
  facet_wrap(~Species) +
  theme_bw() 
#######################################
## Leaf
pnps_cg_l <- read.table("InputData/dNdSpiNpiS_CG_L.txt", header = T, sep = "\t")
head(pnps_cg_l)
pnps_cg_l$CG_L_pNpS <- 0
nna <- pnps_cg_l$CG_piS != 0
pnps_cg_l[nna,]$CG_L_pNpS <- pnps_cg_l[nna,]$CG_piN / pnps_cg_l[nna,]$CG_piS
pnps_cg_l_1 <- pnps_cg_l[, c("Contig_name", "CG_piN", "CG_piS", "CG_L_pNpS")]
row.names(pnps_cg_l_1) <- pnps_cg_l_1$Contig_name
head(pnps_cg_l_1)
# overlap with all genes in RNAseq
all_L <- read.table("InputData/TMM_Diploids_L.txt", header = T, sep = "\t")
row.names(all_L) <- gsub(".g", "", row.names(all_L))
head(all_L)
all_pnps_cg_l <- merge(all_L, pnps_cg_l_1, by = 0) # by row names
head(all_pnps_cg_l)
row.names(all_pnps_cg_l) <- all_pnps_cg_l$Row.names
all_pnps_cg_l$Row.names <- NULL
dim(all_pnps_cg_l)
all_pnps_cg_l <- all_pnps_cg_l[, c("Contig_name", "CG_piN", "CG_piS", "CG_L_pNpS")]
head(all_pnps_cg_l)
# Overlap with coDEGs_CRCO
coDEG_L <- read.table("OutputData/coDEGs_CRCO_L.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_L) <- gsub(".g", "", row.names(coDEG_L))
coDEG_L$coDEGs <- "coDEGs"
head(coDEG_L)
codeg_pnps_cg_l <- merge(all_pnps_cg_l, coDEG_L, by = 0, all = T)
codeg_pnps_cg_l[is.na(codeg_pnps_cg_l$DE_1),]$coDEGs <- "Others"
codeg_pnps_cg_l$Tissue <- "Leaf"
head(codeg_pnps_cg_l)
mean(codeg_pnps_cg_l[codeg_pnps_cg_l$coDEGs == "coDEGs",]$CG_L_pNpS, na.rm = T)
mean(codeg_pnps_cg_l[codeg_pnps_cg_l$coDEGs == "Others",]$CG_L_pNpS, na.rm = T)

# CO_L
pnps_co_l <- read.table("InputData/dNdSpiNpiS_CO_L.txt", header = T, sep = "\t")
head(pnps_co_l)
pnps_co_l$CO_L_pNpS <- 0
nna <- pnps_co_l$CO_piS != 0
pnps_co_l[nna,]$CO_L_pNpS <- pnps_co_l[nna,]$CO_piN / pnps_co_l[nna,]$CO_piS
pnps_co_l_1 <- pnps_co_l[, c("Contig_name", "CO_piN", "CO_piS", "CO_L_pNpS")]
row.names(pnps_co_l_1) <- pnps_co_l_1$Contig_name
head(pnps_co_l_1)
# overlap with all genes in RNAseq
all_pnps_co_l <- merge(all_L, pnps_co_l_1, by = 0) # by row names
head(all_pnps_co_l)
row.names(all_pnps_co_l) <- all_pnps_co_l$Row.names
all_pnps_co_l$Row.names <- NULL
dim(all_pnps_co_l)
all_pnps_co_l <- all_pnps_co_l[, c("Contig_name", "CO_piN", "CO_piS", "CO_L_pNpS")]
head(all_pnps_co_l)
# Overlap with coDEGs_CRCO
coDEG_L <- read.table("OutputData/coDEGs_CRCO_L.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_L) <- gsub(".g", "", row.names(coDEG_L))
coDEG_L$coDEGs <- "coDEGs"
head(coDEG_L)
codeg_pnps_co_l <- merge(all_pnps_co_l, coDEG_L, by = 0, all = T)
codeg_pnps_co_l[is.na(codeg_pnps_co_l$DE_1),]$coDEGs <- "Others"
codeg_pnps_co_l$Tissue <- "Leaf"
head(codeg_pnps_co_l)
mean(codeg_pnps_co_l[codeg_pnps_co_l$coDEGs == "coDEGs",]$CO_L_pNpS, na.rm = T)
mean(codeg_pnps_co_l[codeg_pnps_co_l$coDEGs == "Others",]$CO_L_pNpS, na.rm = T)

### CR
pnps_cr_l <- read.table("InputData/dNdSpiNpiS_CR_L.txt", header = T, sep = "\t")
head(pnps_cr_l)
pnps_cr_l$CR_L_pNpS <- 0
nna <- pnps_cr_l$CR_piS != 0
pnps_cr_l[nna,]$CR_L_pNpS <- pnps_cr_l[nna,]$CR_piN / pnps_cr_l[nna,]$CR_piS
pnps_cr_l_1 <- pnps_cr_l[, c("Contig_name", "CR_piN", "CR_piS", "CR_L_pNpS")]
row.names(pnps_cr_l_1) <- pnps_cr_l_1$Contig_name
head(pnps_cr_l_1)
# overlap with all genes in RNAseq
all_pnps_cr_l <- merge(all_L, pnps_cr_l_1, by = 0) # by row names
head(all_pnps_cr_l)
row.names(all_pnps_cr_l) <- all_pnps_cr_l$Row.names
all_pnps_cr_l$Row.names <- NULL
dim(all_pnps_cr_l)
all_pnps_cr_l <- all_pnps_cr_l[, c("Contig_name", "CR_piN", "CR_piS", "CR_L_pNpS")]
head(all_pnps_cr_l)
# Overlap with coDEGs_CRCO
coDEG_L <- read.table("OutputData/coDEGs_CRCO_L.txt", header = T, sep = "\t", row.names = "Gene")
row.names(coDEG_L) <- gsub(".g", "", row.names(coDEG_L))
coDEG_L$coDEGs <- "coDEGs"
head(coDEG_L)
codeg_pnps_cr_l <- merge(all_pnps_cr_l, coDEG_L, by = 0, all = T)
codeg_pnps_cr_l[is.na(codeg_pnps_cr_l$DE_1),]$coDEGs <- "Others"
codeg_pnps_cr_l$Tissue <- "Leaf"
head(codeg_pnps_cr_l)
median(codeg_pnps_cr_l[codeg_pnps_cr_l$coDEGs == "coDEGs" & codeg_pnps_cr_l$CR_L_pNpS != 0,]$CR_L_pNpS, na.rm = T)
median(codeg_pnps_cr_l[codeg_pnps_cr_l$coDEGs == "Others" & codeg_pnps_cr_l$CR_L_pNpS != 0,]$CR_L_pNpS, na.rm = T)
mean(codeg_pnps_cr_l[codeg_pnps_cr_l$coDEGs == "coDEGs" & codeg_pnps_cr_l$CR_L_pNpS != 0,]$CR_L_pNpS, na.rm = T)
mean(codeg_pnps_cr_l[codeg_pnps_cr_l$coDEGs == "Others" & codeg_pnps_cr_l$CR_L_pNpS != 0,]$CR_L_pNpS, na.rm = T)

# combine CG, CO, CR
head(codeg_pnps_cg_l)
cg_pnps_l <- codeg_pnps_cg_l[,c("Contig_name","CG_piN", "CG_piS", "CG_L_pNpS", "Tissue", "coDEGs")]
cg_pnps_l$Species <- "CG"
cg_pnps_l <- cg_pnps_l[!(is.na(cg_pnps_l$CG_L_pNpS)),]
names(cg_pnps_l) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(cg_pnps_l)

co_pnps_l <- codeg_pnps_co_l[,c("Contig_name","CO_piN", "CO_piS", "CO_L_pNpS", "Tissue", "coDEGs")]
co_pnps_l$Species <- "CO"
co_pnps_l <- co_pnps_l[!(is.na(co_pnps_l$CO_L_pNpS)),]
names(co_pnps_l) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(co_pnps_l)

cr_pnps_l <- codeg_pnps_cr_l[,c("Contig_name","CR_piN", "CR_piS", "CR_L_pNpS", "Tissue", "coDEGs")]
cr_pnps_l$Species <- "CR"
cr_pnps_l <- cr_pnps_l[!(is.na(cr_pnps_l$CR_L_pNpS)),]
names(cr_pnps_l) <- c("ID", "piN", "piS", "piNpiS", "Tissue", "coDEGs", "Species")
head(cr_pnps_l)

pnps_l <- rbind(cg_pnps_l, co_pnps_l, cr_pnps_l)
head(pnps_l)

my_comparisons <- list( c("coDEGs", "Others"))

ggplot(pnps_l, aes(x=coDEGs, y=piNpiS, color=coDEGs)) + 
  geom_boxplot() +
  scale_color_manual(values=c("blue","orange")) +
  facet_grid( ~ Species) +
  stat_compare_means(label.y = 2.5) +
  stat_summary(fun = mean, geom = "point", shape=5, size=1.5, color = "black") +
  #stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
  ylim(0,2.5) +
  theme_bw() +
  labs(x = "Leaves", y = "piN / piS")

ggplot(pnps_l[pnps_l$piNpiS != 0,], aes(x=piNpiS, fill=coDEGs)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_lill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +
  xlim(0,2.5) +
  facet_wrap(~Species) +
  theme_bw() 
