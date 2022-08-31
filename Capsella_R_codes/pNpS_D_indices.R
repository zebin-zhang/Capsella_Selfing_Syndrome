###############################################################
#        piNpiS       Dco_Dcr
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
for (tissue in c("F", "L", "R")) { # Loop_1_tissue_start
  # load all genes in RNAseq
  all_gene <- read.table(paste("InputData/TMM_Diploids_", tissue, ".txt", sep = ""), header = T, sep = "\t")
  row.names(all_gene) <- gsub(".g", "", row.names(all_gene))
  CGColumns <- str_detect(names(all_gene), "CG") # detect which columns contain the "CG" 
  CRColumns <- str_detect(names(all_gene), "CR") # detect which columns contain the "CR" 
  COColumns <- str_detect(names(all_gene), "CO") # detect which columns contain the "CO" 
  
  all_gene$sdCG <- rowSds(as.matrix(all_gene[, CGColumns], na.rm = T))
  all_gene$sdCR <- rowSds(as.matrix(all_gene[, CRColumns], na.rm = T))
  all_gene$sdCO <- rowSds(as.matrix(all_gene[, COColumns], na.rm = T))
  
  all_gene$MeanCG <- rowMeans(all_gene[, CGColumns], na.rm = T)
  all_gene$MeanCR <- rowMeans(all_gene[, CRColumns], na.rm = T)
  all_gene$MeanCO <- rowMeans(all_gene[, COColumns], na.rm = T)
  
  all_gene$cvCG <- all_gene$sdCG / all_gene$MeanCG * 100
  all_gene$cvCR <- all_gene$sdCR / all_gene$MeanCR * 100
  all_gene$cvCO <- all_gene$sdCO / all_gene$MeanCO * 100
  head(all_gene)
  dim(all_gene)
  # load D indices
  # D_CR & D_CO
  if (tissue == "F") {
    Dindex <- read.table("OutputData/Flower_D_indcies_all_genes.txt",header = T, sep = "\t", row.names = 1)
  }else if(tissue == "L"){
    Dindex <- read.table("OutputData/Leaf_D_indcies_all_genes.txt",header = T, sep = "\t", row.names = 1)
  }else{
    Dindex <- read.table("OutputData/Root_D_indcies_all_genes.txt",header = T, sep = "\t", row.names = 1)
  }
  row.names(Dindex) <- gsub(".g", "", row.names(Dindex))
  head(Dindex)
  # merge all_gene with D indices
  all <- merge(all_gene, Dindex, by=0)
  row.names(all) <- all$Row.names
  all$Row.names <- NULL
  head(all)
  dim(all)
  # load pnps files
  for (species in c("CG", "CR", "CO")) {# Loop_2_start
    pnps <- read.table(paste("InputData/dNdSpiNpiS_", species, "_", tissue, ".txt", sep = ""), header = T, sep = "\t", row.names = 1)
    if (species == "CG") {# if_1 start
      nna <- pnps$CG_piS != 0
      pnps <- pnps[nna,] # remove genes with piS = 0.
      pnps_1 <- pnps[, c("CG_nb_complete_site", "CG_piN", "CG_piS")]
    }else if (species == "CR") {
      nna <- pnps$CR_piS != 0
      pnps <- pnps[nna,] # remove genes with piS = 0.
      pnps_1 <- pnps[, c("CR_nb_complete_site", "CR_piN", "CR_piS")]
    }else {
      nna <- pnps$CO_piS != 0
      pnps <- pnps[nna,] # remove genes with piS = 0.
      pnps_1 <- pnps[, c("CO_nb_complete_site", "CO_piN", "CO_piS")]
    }# if_1 end
    
    names(pnps_1) <- c("nb_complete_site", "piN", "piS") # rename pnps1
    # head(pnps)
    # head(pnps_1)
    # dim(pnps_1)
    # merge all_D_indices with pnps
    all_pnps <- merge(all, pnps_1, by = 0) # by row names
    row.names(all_pnps) <- all_pnps$Row.names
    all_pnps$Row.names <- NULL
    head(all_pnps)
    dim(all_pnps)
    
    # all_pnps <- all_pnps[,c("Tissue", "MeanCG", "MeanCR", "MeanCO", "nb_complete_site", "piN", "piS", "D_CR", "D_CO")]
    all_pnps$Species <- species
    
    # if (tissue == "F") {
    #   all_pnps$Tissue <- "Flower"
    # }else if(tissue == "L"){
    #   all_pnps$Tissue <- "Leaf"
    # }else{
    #   all_pnps$Tissue <- "Root"
    # }
    
    all_pnps$ncs_piN <- all_pnps$nb_complete_site * all_pnps$piN # calculate ncs_piN, ncs_piS
    all_pnps$ncs_piS <- all_pnps$nb_complete_site * all_pnps$piS
    all_pnps <- all_pnps[all_pnps$ncs_piN >= 0,]
    assign(paste0("all_pnps_", species, "_", tissue ), all_pnps)
    head(all_pnps)
    write.table(all_pnps, file = paste0("pNpS_data/all_pnps_", species, "_", tissue, ".txt"), append = F, quote = F, sep = "\t", row.names = T)
  }# Loop_2_end
  
}# Loop_1_tissue_end


# When you consider each gene separately you have a huge variance (because most genes have only a few SNPs).
# You compute the ratio mean(piN) / mean(piS) for the categories you want to compare:
# pN_pS = Sum (nb_complete_site x piN) / Sum (nb_complete_site x piS)
head(all_pnps_CG_F)
dim(all_pnps_CG_F)

# pNpS in all genes
value_pnps_CG_F_all <- sum(all_pnps_CG_F$ncs_piN) / sum(all_pnps_CG_F$ncs_piS )
value_pnps_CR_F_all <- sum(all_pnps_CR_F$ncs_piN) / sum(all_pnps_CR_F$ncs_piS )
value_pnps_CO_F_all <- sum(all_pnps_CO_F$ncs_piN) / sum(all_pnps_CO_F$ncs_piS )

value_pnps_CG_L_all <- sum(all_pnps_CG_L$ncs_piN) / sum(all_pnps_CG_L$ncs_piS )
value_pnps_CR_L_all <- sum(all_pnps_CR_L$ncs_piN) / sum(all_pnps_CR_L$ncs_piS )
value_pnps_CO_L_all <- sum(all_pnps_CO_L$ncs_piN) / sum(all_pnps_CO_L$ncs_piS )

value_pnps_CG_R_all <- sum(all_pnps_CG_R$ncs_piN) / sum(all_pnps_CG_R$ncs_piS )
value_pnps_CR_R_all <- sum(all_pnps_CR_R$ncs_piN) / sum(all_pnps_CR_R$ncs_piS )
value_pnps_CO_R_all <- sum(all_pnps_CO_R$ncs_piN) / sum(all_pnps_CO_R$ncs_piS )

# CG_F
head(all_pnps_CG_F)
neg_all <- all_pnps_CG_F$D_CR < 0   & all_pnps_CG_F$D_CO < 0
pos_all <- all_pnps_CG_F$D_CR > 0   & all_pnps_CG_F$D_CO > 0
sum(all_pnps_CG_F[neg_all,]$ncs_piN) / sum(all_pnps_CG_F[neg_all,]$ncs_piS)
sum(all_pnps_CG_F[pos_all,]$ncs_piN) / sum(all_pnps_CG_F[pos_all,]$ncs_piS)
# CR_F
head(all_pnps_CR_F)
neg_all <- all_pnps_CR_F$D_CR < 0   & all_pnps_CR_F$D_CO < 0
pos_all <- all_pnps_CR_F$D_CR > 0   & all_pnps_CR_F$D_CO > 0
sum(all_pnps_CR_F[neg_all,]$ncs_piN) / sum(all_pnps_CR_F[neg_all,]$ncs_piS)
sum(all_pnps_CR_F[pos_all,]$ncs_piN) / sum(all_pnps_CR_F[pos_all,]$ncs_piS)
# CO_F
head(all_pnps_CO_F)
neg_all <- all_pnps_CO_F$D_CR < 0   & all_pnps_CO_F$D_CO < 0
pos_all <- all_pnps_CO_F$D_CR > 0   & all_pnps_CO_F$D_CO > 0
sum(all_pnps_CO_F[neg_all,]$ncs_piN) / sum(all_pnps_CO_F[neg_all,]$ncs_piS)
sum(all_pnps_CO_F[pos_all,]$ncs_piN) / sum(all_pnps_CO_F[pos_all,]$ncs_piS)
# CG_L
head(all_pnps_CG_L)
neg_all <- all_pnps_CG_L$D_CR < 0   & all_pnps_CG_L$D_CO < 0
pos_all <- all_pnps_CG_L$D_CR > 0   & all_pnps_CG_L$D_CO > 0
sum(all_pnps_CG_L[neg_all,]$ncs_piN) / sum(all_pnps_CG_L[neg_all,]$ncs_piS)
sum(all_pnps_CG_L[pos_all,]$ncs_piN) / sum(all_pnps_CG_L[pos_all,]$ncs_piS)
# CR_L
head(all_pnps_CR_L)
neg_all <- all_pnps_CR_L$D_CR < 0   & all_pnps_CR_L$D_CO < 0
pos_all <- all_pnps_CR_L$D_CR > 0   & all_pnps_CR_L$D_CO > 0
sum(all_pnps_CR_L[neg_all,]$ncs_piN) / sum(all_pnps_CR_L[neg_all,]$ncs_piS)
sum(all_pnps_CR_L[pos_all,]$ncs_piN) / sum(all_pnps_CR_L[pos_all,]$ncs_piS)
# CO_L
head(all_pnps_CO_L)
neg_all <- all_pnps_CO_L$D_CR < 0   & all_pnps_CO_L$D_CO < 0
pos_all <- all_pnps_CO_L$D_CR > 0   & all_pnps_CO_L$D_CO > 0
sum(all_pnps_CO_L[neg_all,]$ncs_piN) / sum(all_pnps_CO_L[neg_all,]$ncs_piS)
sum(all_pnps_CO_L[pos_all,]$ncs_piN) / sum(all_pnps_CO_L[pos_all,]$ncs_piS)
# CG_R
head(all_pnps_CG_R)
neg_all <- all_pnps_CG_R$D_CR < 0   & all_pnps_CG_R$D_CO < 0
pos_all <- all_pnps_CG_R$D_CR > 0   & all_pnps_CG_R$D_CO > 0
sum(all_pnps_CG_R[neg_all,]$ncs_piN) / sum(all_pnps_CG_R[neg_all,]$ncs_piS)
sum(all_pnps_CG_R[pos_all,]$ncs_piN) / sum(all_pnps_CG_R[pos_all,]$ncs_piS)
# CR_R
head(all_pnps_CR_R)
neg_all <- all_pnps_CR_R$D_CR < 0   & all_pnps_CR_R$D_CO < 0
pos_all <- all_pnps_CR_R$D_CR > 0   & all_pnps_CR_R$D_CO > 0
sum(all_pnps_CR_R[neg_all,]$ncs_piN) / sum(all_pnps_CR_R[neg_all,]$ncs_piS)
sum(all_pnps_CR_R[pos_all,]$ncs_piN) / sum(all_pnps_CR_R[pos_all,]$ncs_piS)
# CO_R
head(all_pnps_CO_R)
neg_all <- all_pnps_CO_R$D_CR < 0   & all_pnps_CO_R$D_CO < 0
pos_all <- all_pnps_CO_R$D_CR > 0   & all_pnps_CO_R$D_CO > 0
sum(all_pnps_CO_R[neg_all,]$ncs_piN) / sum(all_pnps_CO_R[neg_all,]$ncs_piS)
sum(all_pnps_CO_R[pos_all,]$ncs_piN) / sum(all_pnps_CO_R[pos_all,]$ncs_piS)
# pNpS in all Neg- and Pos-overlap of Dco Dcr
for (tissue in c("F", "L", "R")) { # Loop_1_tissue_start
  for (species in c("CG", "CR", "CO")) {# Loop_2_start
    all_pnps <- read.table(paste0("pNpS_data/all_pnps_", species, "_", tissue, ".txt"), header = T, sep = "\t")
    head(all_pnps)
    for (drange in c(0.1, 0.3, 0.5, 0.7, 0.9)){#loop_3_Start
      neg_d <- all_pnps$D_CO < -drange & all_pnps$D_CR < -drange
      pos_d <- all_pnps$D_CO >  drange & all_pnps$D_CR >  drange
      pnps_neg_d <- sum(all_pnps[neg_d,]$ncs_piN) / sum(all_pnps[neg_d,]$ncs_piS)
      pnps_pos_d <- sum(all_pnps[pos_d,]$ncs_piN) / sum(all_pnps[pos_d,]$ncs_piS)
      
      assign(paste0("neg_", drange, "_pnps_", species, "_", tissue), pnps_neg_d)
      assign(paste0("pos_", drange, "_pnps_", species, "_", tissue), pnps_pos_d)
    }#loop_3_end
  }#loop_2_end
} # loop_1_end

# CG_F
head(all_pnps_CG_F)
neg_all <- all_pnps_CG_F$D_CR < 0 & all_pnps_CG_F$D_CO < 0
neg_0.1 <- all_pnps_CG_F$D_CR < -0.1 & all_pnps_CG_F$D_CO < -0.1
neg_0.3 <- all_pnps_CG_F$D_CR < -0.3 & all_pnps_CG_F$D_CO < -0.3
neg_0.5 <- all_pnps_CG_F$D_CR < -0.5 & all_pnps_CG_F$D_CO < -0.5
neg_0.7 <- all_pnps_CG_F$D_CR < -0.7 & all_pnps_CG_F$D_CO < -0.7
neg_0.9 <- all_pnps_CG_F$D_CR < -0.9 & all_pnps_CG_F$D_CO < -0.9

pos_all <- all_pnps_CG_F$D_CR > 0   & all_pnps_CG_F$D_CO > 0
pos_0.1 <- all_pnps_CG_F$D_CR > 0.1 & all_pnps_CG_F$D_CO > 0.1
pos_0.3 <- all_pnps_CG_F$D_CR > 0.3 & all_pnps_CG_F$D_CO > 0.3
pos_0.5 <- all_pnps_CG_F$D_CR > 0.5 & all_pnps_CG_F$D_CO > 0.5
pos_0.7 <- all_pnps_CG_F$D_CR > 0.7 & all_pnps_CG_F$D_CO > 0.7
pos_0.9 <- all_pnps_CG_F$D_CR > 0.9 & all_pnps_CG_F$D_CO > 0.9

head(all_pnps_CG_F[neg_all,])
head(all_pnps_CG_F[neg_0.9,])

sum(all_pnps_CG_F[neg_all,]$ncs_piN) / sum(all_pnps_CG_F[neg_all,]$ncs_piS)
sum(all_pnps_CG_F[neg_0.1,]$ncs_piN) / sum(all_pnps_CG_F[neg_0.1,]$ncs_piS)
sum(all_pnps_CG_F[neg_0.3,]$ncs_piN) / sum(all_pnps_CG_F[neg_0.3,]$ncs_piS)
sum(all_pnps_CG_F[neg_0.5,]$ncs_piN) / sum(all_pnps_CG_F[neg_0.5,]$ncs_piS)
sum(all_pnps_CG_F[neg_0.7,]$ncs_piN) / sum(all_pnps_CG_F[neg_0.7,]$ncs_piS)
sum(all_pnps_CG_F[neg_0.9,]$ncs_piN) / sum(all_pnps_CG_F[neg_0.9,]$ncs_piS)

sum(all_pnps_CG_F[pos_all,]$ncs_piN) / sum(all_pnps_CG_F[pos_all,]$ncs_piS)
sum(all_pnps_CG_F[pos_0.1,]$ncs_piN) / sum(all_pnps_CG_F[pos_0.1,]$ncs_piS)
sum(all_pnps_CG_F[pos_0.3,]$ncs_piN) / sum(all_pnps_CG_F[pos_0.3,]$ncs_piS)
sum(all_pnps_CG_F[pos_0.5,]$ncs_piN) / sum(all_pnps_CG_F[pos_0.5,]$ncs_piS)
sum(all_pnps_CG_F[pos_0.7,]$ncs_piN) / sum(all_pnps_CG_F[pos_0.7,]$ncs_piS)
sum(all_pnps_CG_F[pos_0.9,]$ncs_piN) / sum(all_pnps_CG_F[pos_0.9,]$ncs_piS)

# CR_F
head(all_pnps_CR_F)
dim(all_pnps_CR_F)
neg_all <- all_pnps_CR_F$D_CR < 0 & all_pnps_CR_F$D_CO < 0
neg_0.1 <- all_pnps_CR_F$D_CR < -0.1 & all_pnps_CR_F$D_CO < -0.1
neg_0.3 <- all_pnps_CR_F$D_CR < -0.3 & all_pnps_CR_F$D_CO < -0.3
neg_0.5 <- all_pnps_CR_F$D_CR < -0.5 & all_pnps_CR_F$D_CO < -0.5
neg_0.7 <- all_pnps_CR_F$D_CR < -0.7 & all_pnps_CR_F$D_CO < -0.7
neg_0.9 <- all_pnps_CR_F$D_CR < -0.9 & all_pnps_CR_F$D_CO < -0.9

pos_all <- all_pnps_CR_F$D_CR > 0   & all_pnps_CR_F$D_CO > 0
pos_0.1 <- all_pnps_CR_F$D_CR > 0.1 & all_pnps_CR_F$D_CO > 0.1
pos_0.3 <- all_pnps_CR_F$D_CR > 0.3 & all_pnps_CR_F$D_CO > 0.3
pos_0.5 <- all_pnps_CR_F$D_CR > 0.5 & all_pnps_CR_F$D_CO > 0.5
pos_0.7 <- all_pnps_CR_F$D_CR > 0.7 & all_pnps_CR_F$D_CO > 0.7
pos_0.9 <- all_pnps_CR_F$D_CR > 0.9 & all_pnps_CR_F$D_CO > 0.9

head(all_pnps_CR_F[neg_all,])
head(all_pnps_CR_F[neg_0.9,])

sum(all_pnps_CR_F[neg_all,]$ncs_piN) / sum(all_pnps_CR_F[neg_all,]$ncs_piS)
sum(all_pnps_CR_F[neg_0.1,]$ncs_piN) / sum(all_pnps_CR_F[neg_0.1,]$ncs_piS)
sum(all_pnps_CR_F[neg_0.3,]$ncs_piN) / sum(all_pnps_CR_F[neg_0.3,]$ncs_piS)
sum(all_pnps_CR_F[neg_0.5,]$ncs_piN) / sum(all_pnps_CR_F[neg_0.5,]$ncs_piS)
sum(all_pnps_CR_F[neg_0.7,]$ncs_piN) / sum(all_pnps_CR_F[neg_0.7,]$ncs_piS)
sum(all_pnps_CR_F[neg_0.9,]$ncs_piN) / sum(all_pnps_CR_F[neg_0.9,]$ncs_piS)

sum(all_pnps_CR_F[pos_all,]$ncs_piN) / sum(all_pnps_CR_F[pos_all,]$ncs_piS)
sum(all_pnps_CR_F[pos_0.1,]$ncs_piN) / sum(all_pnps_CR_F[pos_0.1,]$ncs_piS)
sum(all_pnps_CR_F[pos_0.3,]$ncs_piN) / sum(all_pnps_CR_F[pos_0.3,]$ncs_piS)
sum(all_pnps_CR_F[pos_0.5,]$ncs_piN) / sum(all_pnps_CR_F[pos_0.5,]$ncs_piS)
sum(all_pnps_CR_F[pos_0.7,]$ncs_piN) / sum(all_pnps_CR_F[pos_0.7,]$ncs_piS)
sum(all_pnps_CR_F[pos_0.9,]$ncs_piN) / sum(all_pnps_CR_F[pos_0.9,]$ncs_piS)


# bootstrap CI 95%
pnps <- function(d, i){
  d2 <- d[i,]
  pn <- d2$nb_complete_site * d2$piN
  ps <- d2$nb_complete_site * d2$piS
  return(sum(pn) / sum(ps))
}

bootstrap_all_pnps_cg_f <- boot(all_pnps_CG_F, pnps, R=1000)
summary(bootstrap_all_pnps_cg_f)
# boot.ci(boot.out=bootstrap_all_pnps_cg_f,type=c("norm","basic","perc","bca"))
sort(bootstrap_all_pnps_cg_f$t)[c(25,975)]

bootstrap_all_pnps_cr_f <- boot(all_pnps_CR_F, pnps, R=1000)
summary(bootstrap_all_pnps_cr_f)
sort(bootstrap_all_pnps_cr_f$t)[c(25,975)]
# boot.ci(boot.out=bootstrap_all_pnps_cr_f,type=c("norm","basic","perc","bca"))

bootstrap_all_pnps_co_f <- boot(all_pnps_CO_F, pnps, R=1000)
summary(bootstrap_all_pnps_co_f)
sort(bootstrap_all_pnps_co_f$t)[c(25,975)]
# boot.ci(boot.out=bootstrap_all_pnps_co_f,type=c("norm","basic","perc","bca"))


































