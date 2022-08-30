###############################################################
#        coefficient of variation       
###############################################################
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
  assign(paste0("all_CV_Dindices_", tissue ), all)
  write.table(all, file = paste0("coefficient_of_variation/all_CV_Dindices_",  tissue, ".txt"), append = F, quote = F, sep = "\t", row.names = T)
}# Loop_1_tissue_end

# all genes
mean(all_CV_Dindices_F$cvCG, na.rm = T)
mean(all_CV_Dindices_F$cvCR, na.rm = T)
mean(all_CV_Dindices_F$cvCO, na.rm = T)

mean(all_CV_Dindices_L$cvCG, na.rm = T)
mean(all_CV_Dindices_L$cvCR, na.rm = T)
mean(all_CV_Dindices_L$cvCO, na.rm = T)

mean(all_CV_Dindices_R$cvCG, na.rm = T)
mean(all_CV_Dindices_R$cvCR, na.rm = T)
mean(all_CV_Dindices_R$cvCO, na.rm = T)

##################### F - Dco 
neg_all <- all_CV_Dindices_F$D_CO < 0
neg_0.1 <- all_CV_Dindices_F$D_CO < -0.1
neg_0.3 <- all_CV_Dindices_F$D_CO < -0.3
neg_0.5 <- all_CV_Dindices_F$D_CO < -0.5
neg_0.7 <- all_CV_Dindices_F$D_CO < -0.7
neg_0.9 <- all_CV_Dindices_F$D_CO < -0.9

pos_all <- all_CV_Dindices_F$D_CO > 0
pos_0.1 <- all_CV_Dindices_F$D_CO > 0.1
pos_0.3 <- all_CV_Dindices_F$D_CO > 0.3
pos_0.5 <- all_CV_Dindices_F$D_CO > 0.5
pos_0.7 <- all_CV_Dindices_F$D_CO > 0.7
pos_0.9 <- all_CV_Dindices_F$D_CO > 0.9
# CG NEG
mean(all_CV_Dindices_F[neg_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[neg_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[neg_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[neg_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[neg_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[neg_0.9,]$cvCG, na.rm = T)
# CG POS
mean(all_CV_Dindices_F[pos_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[pos_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[pos_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[pos_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[pos_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_F[pos_0.9,]$cvCG, na.rm = T)
# CR NEG
mean(all_CV_Dindices_F[neg_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[neg_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[neg_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[neg_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[neg_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[neg_0.9,]$cvCR, na.rm = T)
# CR POS
mean(all_CV_Dindices_F[pos_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[pos_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[pos_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[pos_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[pos_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_F[pos_0.9,]$cvCR, na.rm = T)
# CO NEG
mean(all_CV_Dindices_F[neg_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[neg_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[neg_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[neg_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[neg_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[neg_0.9,]$cvCO, na.rm = T)
# CO POS
mean(all_CV_Dindices_F[pos_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[pos_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[pos_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[pos_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[pos_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_F[pos_0.9,]$cvCO, na.rm = T)

##################### L - Dco  overlaps
neg_all <-  all_CV_Dindices_L$D_CO < 0
neg_0.1 <-  all_CV_Dindices_L$D_CO < -0.1
neg_0.3 <-  all_CV_Dindices_L$D_CO < -0.3
neg_0.5 <-  all_CV_Dindices_L$D_CO < -0.5
neg_0.7 <-  all_CV_Dindices_L$D_CO < -0.7
neg_0.9 <-  all_CV_Dindices_L$D_CO < -0.9

pos_all <- all_CV_Dindices_L$D_CO > 0
pos_0.1 <- all_CV_Dindices_L$D_CO > 0.1
pos_0.3 <- all_CV_Dindices_L$D_CO > 0.3
pos_0.5 <- all_CV_Dindices_L$D_CO > 0.5
pos_0.7 <- all_CV_Dindices_L$D_CO > 0.7
pos_0.9 <- all_CV_Dindices_L$D_CO > 0.9
# CG NEG
mean(all_CV_Dindices_L[neg_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[neg_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[neg_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[neg_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[neg_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[neg_0.9,]$cvCG, na.rm = T)
# CG POS
mean(all_CV_Dindices_L[pos_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[pos_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[pos_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[pos_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[pos_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_L[pos_0.9,]$cvCG, na.rm = T)
# CR NEG
mean(all_CV_Dindices_L[neg_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[neg_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[neg_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[neg_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[neg_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[neg_0.9,]$cvCR, na.rm = T)
# CR POS
mean(all_CV_Dindices_L[pos_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[pos_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[pos_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[pos_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[pos_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_L[pos_0.9,]$cvCR, na.rm = T)
# CO NEG
mean(all_CV_Dindices_L[neg_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[neg_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[neg_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[neg_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[neg_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[neg_0.9,]$cvCO, na.rm = T)
# CO POS
mean(all_CV_Dindices_L[pos_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[pos_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[pos_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[pos_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[pos_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_L[pos_0.9,]$cvCO, na.rm = T)

##################### R - Dco Dcr overlaps
neg_all <-  all_CV_Dindices_R$D_CO < 0
neg_0.1 <-  all_CV_Dindices_R$D_CO < -0.1
neg_0.3 <-  all_CV_Dindices_R$D_CO < -0.3
neg_0.5 <-  all_CV_Dindices_R$D_CO < -0.5
neg_0.7 <-  all_CV_Dindices_R$D_CO < -0.7
neg_0.9 <-  all_CV_Dindices_R$D_CO < -0.9

pos_all <- all_CV_Dindices_R$D_CO > 0
pos_0.1 <- all_CV_Dindices_R$D_CO > 0.1
pos_0.3 <- all_CV_Dindices_R$D_CO > 0.3
pos_0.5 <- all_CV_Dindices_R$D_CO > 0.5
pos_0.7 <- all_CV_Dindices_R$D_CO > 0.7
pos_0.9 <- all_CV_Dindices_R$D_CO > 0.9
# CG NEG
mean(all_CV_Dindices_R[neg_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[neg_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[neg_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[neg_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[neg_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[neg_0.9,]$cvCG, na.rm = T)
# CG POS
mean(all_CV_Dindices_R[pos_all,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[pos_0.1,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[pos_0.3,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[pos_0.5,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[pos_0.7,]$cvCG, na.rm = T)
mean(all_CV_Dindices_R[pos_0.9,]$cvCG, na.rm = T)
# CR NEG
mean(all_CV_Dindices_R[neg_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[neg_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[neg_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[neg_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[neg_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[neg_0.9,]$cvCR, na.rm = T)
# CR POS
mean(all_CV_Dindices_R[pos_all,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[pos_0.1,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[pos_0.3,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[pos_0.5,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[pos_0.7,]$cvCR, na.rm = T)
mean(all_CV_Dindices_R[pos_0.9,]$cvCR, na.rm = T)
# CO NEG
mean(all_CV_Dindices_R[neg_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[neg_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[neg_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[neg_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[neg_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[neg_0.9,]$cvCO, na.rm = T)
# CO POS
mean(all_CV_Dindices_R[pos_all,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[pos_0.1,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[pos_0.3,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[pos_0.5,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[pos_0.7,]$cvCO, na.rm = T)
mean(all_CV_Dindices_R[pos_0.9,]$cvCO, na.rm = T)











# for (tissue in c("F", "L", "R")) {
#   Reads <- read.table(paste("InputData/TMM_Diploids_", tissue, ".txt", sep = ""), header = T, sep = "\t")
#   head(Reads)
#   CGColumns <- str_detect(names(Reads), "CG") # detect which columns contain the "CG" 
#   CRColumns <- str_detect(names(Reads), "CR") # detect which columns contain the "CR" 
#   COColumns <- str_detect(names(Reads), "CO") # detect which columns contain the "CO" 
#   
#   Reads$MeanCG <- rowMeans(Reads[, CGColumns], na.rm = T)
#   Reads$MeanCR <- rowMeans(Reads[, CRColumns], na.rm = T)
#   Reads$MeanCO <- rowMeans(Reads[, COColumns], na.rm = T)
#   Reads$MeanSelfers <- rowMeans(Reads[, CRColumns | COColumns], na.rm = T)
#   
#   Reads$sdCG <- rowSds(as.matrix(Reads[, CGColumns], na.rm = T))
#   Reads$sdCR <- rowSds(as.matrix(Reads[, CRColumns], na.rm = T))
#   Reads$sdCO <- rowSds(as.matrix(Reads[, COColumns], na.rm = T))
#   Reads$sdSelfers <- rowSds(as.matrix(Reads[, CRColumns | COColumns], na.rm = T))
#   
#   Reads$cvCG <- Reads$sdCG / Reads$MeanCG * 100
#   Reads$cvCR <- Reads$sdCR / Reads$MeanCR * 100
#   Reads$cvCO <- Reads$sdCO / Reads$MeanCO * 100
#   Reads$cvSelfers <- Reads$sdSelfers / Reads$MeanSelfers * 100
#   
#   assign(paste0("Reads_", tissue), Reads)
# }
# head(Reads)
# 
# # Overlap with coDEGs_CRCO
# coDEG_F <- read.table("OutputData/coDEGs_CRCO_F.txt", header = T, sep = "\t", row.names = "Gene")
# head(coDEG_F)
# head(Reads_F)
# Reads_coDEG_F <- merge(Reads_F, coDEG_F, by = 0, all = T)
# head(Reads_coDEG_F)
# Reads_coDEG_F[is.na(Reads_coDEG_F$DE_1),]$DE_1 <- "Others"
# 
# data_F <- Reads_coDEG_F[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1")]
# melt_data_F <- melt(data_F, id = c("DE_1"))
# names(melt_data_F) <- c("coDEGs", "groups", "cv")
# head(melt_data_F)
# melt_data_F$tissue <- "Flower"
# melt_data_F$coDEGs <- factor(melt_data_F$coDEGs, levels = c("Down", "Up", "Others"))
# 
# ggplot(melt_data_F, aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("coDEGs") +
#   ylim(0,300)
#   
# coDEG_L <- read.table("OutputData/coDEGs_CRCO_L.txt", header = T, sep = "\t", row.names = "Gene")
# Reads_coDEG_L <- merge(Reads_L, coDEG_L, by = 0, all = T)
# Reads_coDEG_L[is.na(Reads_coDEG_L$DE_1),]$DE_1 <- "Others"
# data_L <- Reads_coDEG_L[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1")]
# melt_data_L <- melt(data_L, id = c("DE_1"))
# names(melt_data_L) <- c("coDEGs", "groups", "cv")
# head(melt_data_L)
# melt_data_L$tissue <- "Leaf"
# melt_data_L$coDEGs <- factor(melt_data_L$coDEGs, levels = c("Down", "Up", "Others"))
# ggplot(melt_data_L, aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("coDEGs") +
#   ylim(0,300)
# 
# 
# coDEG_R <- read.table("OutputData/coDEGs_CRCO_R.txt", header = T, sep = "\t", row.names = "Gene")
# Reads_coDEG_R <- merge(Reads_R, coDEG_R, by = 0, all = T)
# Reads_coDEG_R[is.na(Reads_coDEG_R$DE_1),]$DE_1 <- "Others"
# data_R <- Reads_coDEG_R[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1")]
# melt_data_R <- melt(data_R, id = c("DE_1"))
# names(melt_data_R) <- c("coDEGs", "groups", "cv")
# head(melt_data_R)
# melt_data_R$tissue <- "Root"
# melt_data_R$coDEGs <- factor(melt_data_R$coDEGs, levels = c("Down", "Up", "Others"))
# 
# my_comparisons <- list( c("cvCG", "cvCR"), c("cvCG", "cvCO"))
# ggplot(melt_data_R[melt_data_R$cv < 300,], aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("coDEGs") +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)
# 
# data_FLR <- rbind(melt_data_F, melt_data_L, melt_data_R)
# head(data_FLR)
# str(data_FLR)
# data_FLR$coDEGs <- recode_factor(data_FLR$coDEGs, Down = "CRCO coDEGs Down", 
#                                 Up = "CRCO coDEGs Up")
# 
# data_FLR$groups <- recode_factor(data_FLR$groups, cvCG = "CG", 
#                                  cvCR = "CR", 
#                                  cvCO = "CO",
#                                  cvSelfers = "Selfers")
# 
# my_comparisons <- list( c("CG", "CR"), c("CG", "CO"))
# 
# ggplot(data_FLR[data_FLR$cv < 300,], aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#     facet_grid(tissue ~ coDEGs) +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
#   ylim(0,375) +
#   theme_bw() 
# 
# data_FLR_1 <- data_FLR[data_FLR$groups != "Selfers",]
# head(data_FLR_1)
# names(data_FLR_1) <- c("coDEGs", "species", "cv", "tissue")
# ggplot(data_FLR_1[data_FLR_1$cv < 300,], aes(x=species, y=cv, color = species)) + 
#   geom_boxplot() +
#   scale_color_manual(values=c("red", "purple", "blue")) +
#   facet_grid(tissue ~ coDEGs) +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
#   ylim(0,375) +
#   theme_bw() +
#   labs(y = "coefficient of variation")
# 
# # overlap with MST genes
# 
# MST <- read.table("OutputData/482_MST_genes_overlap_Slotte_Wozniak.txt", header = T, sep = "\t")
# MST$gene_tpye <- "MST"
# MST$DE_1 <- MST$DE # create new column for CR/CG
# MST[MST$DE == "Up",]$DE_1 <- "Down" # convert CG/CR Up to CR/CG Down
# MST[MST$DE == "Down",]$DE_1 <- "Up" # convert CG/CR Down to CR/CG Up
# MST <- MST[, c("DE_1", "gene_tpye")]
# # row.names(MST) <- gsub(".g", "", row.names(MST))
# head(MST)
# # flower
# Reads_MST_F <- merge(Reads_F, MST, by = 0, all = T)
# head(Reads_MST_F)
# Reads_MST_F[is.na(Reads_MST_F$DE_1),]$DE_1 <- "Others"
# Reads_MST_F[is.na(Reads_MST_F$gene_tpye),]$gene_tpye <- "Others"
# 
# data_MST_F <- Reads_MST_F[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1","gene_tpye")]
# head(data_MST_F)
# melt_data_MST_F <- melt(data_MST_F, id = c("DE_1", "gene_tpye"))
# names(melt_data_MST_F) <- c("DE", "gene_type", "groups", "cv")
# head(melt_data_MST_F)
# melt_data_MST_F$tissue <- "Flower"
# melt_data_MST_F$gene_type <- factor(melt_data_MST_F$gene_type, levels = c("MST", "Others"))
# melt_data_MST_F$DE <- factor(melt_data_MST_F$DE, levels = c("Up", "Down", "Others"))
# 
# ggplot(melt_data_MST_F, aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# 
# ggplot(melt_data_MST_F, aes(x=groups, y=cv, color=gene_type, fill=DE)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# # leaf
# Reads_MST_L <- merge(Reads_L, MST, by = 0, all = T)
# head(Reads_MST_L)
# Reads_MST_L[is.na(Reads_MST_L$DE_1),]$DE_1 <- "Others"
# Reads_MST_L[is.na(Reads_MST_L$gene_tpye),]$gene_tpye <- "Others"
# 
# data_MST_L <- Reads_MST_L[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1","gene_tpye")]
# head(data_MST_L)
# melt_data_MST_L <- melt(data_MST_L, id = c("DE_1", "gene_tpye"))
# names(melt_data_MST_L) <- c("DE", "gene_type", "groups", "cv")
# head(melt_data_MST_L)
# melt_data_MST_L$tissue <- "Leaf"
# melt_data_MST_L$gene_type <- factor(melt_data_MST_L$gene_type, levels = c("MST", "Others"))
# melt_data_MST_L$DE <- factor(melt_data_MST_L$DE, levels = c("Up", "Down", "Others"))
# 
# ggplot(melt_data_MST_L, aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# 
# ggplot(melt_data_MST_L, aes(x=groups, y=cv, color=gene_type, fill=DE)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# # Root
# Reads_MST_R <- merge(Reads_R, MST, by = 0, all = T)
# head(Reads_MST_R)
# Reads_MST_R[is.na(Reads_MST_R$DE_1),]$DE_1 <- "Others"
# Reads_MST_R[is.na(Reads_MST_R$gene_tpye),]$gene_tpye <- "Others"
# 
# data_MST_R <- Reads_MST_R[,c("cvCG","cvCR","cvCO","cvSelfers","DE_1","gene_tpye")]
# head(data_MST_R)
# melt_data_MST_R <- melt(data_MST_R, id = c("DE_1", "gene_tpye"))
# names(melt_data_MST_R) <- c("DE", "gene_type", "groups", "cv")
# head(melt_data_MST_R)
# melt_data_MST_R$tissue <- "Root"
# melt_data_MST_R$gene_type <- factor(melt_data_MST_R$gene_type, levels = c("MST", "Others"))
# melt_data_MST_R$DE <- factor(melt_data_MST_R$DE, levels = c("Up", "Down", "Others"))
# 
# ggplot(melt_data_MST_R, aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# 
# ggplot(melt_data_MST_R, aes(x=groups, y=cv, color=gene_type, fill=DE)) + 
#   geom_boxplot() +
#   facet_wrap("gene_type") +
#   ylim(0,300)
# 
# # Combine FLR
# head(melt_data_MST_F)
# head(melt_data_MST_L)
# head(melt_data_MST_R)
# 
# data_MST_FLR <- rbind(melt_data_MST_F, melt_data_MST_L, melt_data_MST_R)
# head(data_MST_FLR)
# str(data_MST_FLR)
# data_MST_FLR$MST <- data_MST_FLR$DE
# data_MST_FLR$MST <- as.factor(data_MST_FLR$MST, levels=c("Up","Down","Others"))
# 
# data_MST_FLR$MST <- recode_factor(data_MST_FLR$MST, Up = "MST_Up", 
#                                      Down = "MST_Down", 
#                                      Others = "Others")
# 
# data_MST_FLR$groups <- recode_factor(data_MST_FLR$groups, cvCG = "CG", 
#                                  cvCR = "CR", 
#                                  cvCO = "CO",
#                                  cvSelfers = "Selfers")
# 
# my_comparisons <- list( c("CG", "CR"), c("CG", "CO"))
# 
# ggplot(data_MST_FLR[data_MST_FLR$cv < 300,], aes(x=groups, y=cv)) + 
#   geom_boxplot() +
#   facet_grid(tissue ~ MST) +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
#   ylim(0,375) +
#   theme_bw() 
# 
# data_MST_FLR_1 <- data_MST_FLR[data_MST_FLR$groups != "Selfers",]
# head(data_MST_FLR_1)
# names(data_MST_FLR_1) <- c("DE", "gene_type", "species", "cv", "tissue", "MST")
# ggplot(data_MST_FLR_1[data_MST_FLR_1$cv < 300,], aes(x=species, y=cv, color = species)) + 
#   geom_boxplot() +
#   scale_color_manual(values=c("red", "purple", "blue")) +
#   facet_grid(tissue ~ MST) +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
#   ylim(0,375) +
#   theme_bw() +
#   labs(y = "coefficient of variation")
# 
