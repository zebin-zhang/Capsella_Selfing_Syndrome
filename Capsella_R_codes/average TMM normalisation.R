
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

# Install requested libraries from Bioconductor
# try http:// if https:// URLs are not supported
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("edgeR")
# BiocManager::install("limma")
# # Load required libraries 
library(DESeq2)
library(edgeR)
library(limma)
library(Biobase)


#### This file contains reads counts infromation of each gene in all samples
# Phased data
raw.RC <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE_rename.txt", header=T, row.names="gene")
head(raw.RC)
dim(raw.RC)
# CPMnorm <- read.table("Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv", header = T, sep = "\t")
# head(CPMnorm)
# dim(CPMnorm)
##### This file contains information of population and tissue in each sample
Accession <- read.csv("InputData/28genomes_annot_R.txt", header= T, sep = "\t")
head(Accession)
str(Accession)
### Define populations
######################################################
CR_F <- c(Accession$populations=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$populations=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$populations=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$populations=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$populations=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$populations=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$populations=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$populations=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$populations=="CO" & Accession$tissue=="root")

Cbp_Co <- c("ASI_Co", "EUR_Co", "ME_Co", "CASI_Co")
Cbp_Cg <- c("ASI_Cg", "EUR_Cg", "ME_Cg", "CASI_Cg")

Cbp_Co_F <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="flower")
Cbp_Co_L <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="leaf")
Cbp_Co_R <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="root")
Cbp_Cg_F <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="flower")
Cbp_Cg_L <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="leaf")
Cbp_Cg_R <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="root")

ASI_Co_F <- c(Accession$populations=="ASI_Co" & Accession$tissue=="flower")
EUR_Co_F <- c(Accession$populations=="EUR_Co" & Accession$tissue=="flower")
ME_Co_F <- c(Accession$populations=="ME_Co" & Accession$tissue=="flower")
CASI_Co_F <- c(Accession$populations=="CASI_Co" & Accession$tissue=="flower")

ASI_Cg_F <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="flower")
EUR_Cg_F <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="flower")
ME_Cg_F <- c(Accession$populations=="ME_Cg" & Accession$tissue=="flower")
CASI_Cg_F <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="flower")

ASI_Co_L <- c(Accession$populations=="ASI_Co" & Accession$tissue=="leaf")
EUR_Co_L <- c(Accession$populations=="EUR_Co" & Accession$tissue=="leaf")
ME_Co_L <- c(Accession$populations=="ME_Co" & Accession$tissue=="leaf")
CASI_Co_L <- c(Accession$populations=="CASI_Co" & Accession$tissue=="leaf")

ASI_Cg_L <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="leaf")
EUR_Cg_L <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="leaf")
ME_Cg_L <- c(Accession$populations=="ME_Cg" & Accession$tissue=="leaf")
CASI_Cg_L <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="leaf")

ASI_Co_R <- c(Accession$populations=="ASI_Co" & Accession$tissue=="root")
EUR_Co_R <- c(Accession$populations=="EUR_Co" & Accession$tissue=="root")
ME_Co_R <- c(Accession$populations=="ME_Co" & Accession$tissue=="root")
CASI_Co_R <- c(Accession$populations=="CASI_Co" & Accession$tissue=="root")

ASI_Cg_R <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="root")
EUR_Cg_R <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="root")
ME_Cg_R <- c(Accession$populations=="ME_Cg" & Accession$tissue=="root")
CASI_Cg_R <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="root")
#####################################################
# Filter reads
geneMis <- 1 
geneKeep <- cbind((rowSums(raw.RC[,CG_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_F] >= 1, na.rm = T) >= geneMis) &
                    
                    (rowSums(raw.RC[,CG_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_L] >= 1, na.rm = T) >= geneMis) &
                    
                    (rowSums(raw.RC[,CG_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_R] >= 1, na.rm = T) >= geneMis) &
                  ## Run following commands if analysis include Cbp.                     
                    (rowSums(raw.RC[,ASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Cg_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Cg_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
                              
                    (rowSums(raw.RC[,ASI_Co_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Co_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Co_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Co_F] >= 1, na.rm = T) >= geneMis) &
                              
                    (rowSums(raw.RC[,ASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Cg_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Cg_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
                              
                    (rowSums(raw.RC[,ASI_Co_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Co_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Co_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Co_L] >= 1, na.rm = T) >= geneMis) &
                              
                    (rowSums(raw.RC[,ASI_Cg_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Cg_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Cg_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Cg_R] >= 1, na.rm = T) >= geneMis) &
                    
                    (rowSums(raw.RC[,ASI_Co_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,EUR_Co_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,ME_Co_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CASI_Co_R] >= 1, na.rm = T) >= geneMis)
)

raw.RC.filter1 <- raw.RC[geneKeep,]
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.filter1)[1]))

## Filter HSE
# Filter by NA number
# For each gene in each tissue in parental species, at least 5 samples contain none-zero value.
Fkeep <- rowSums(raw.RC.filter1[,c(1:12)] >= 1, na.rm = T) >= 1
Lkeep <- rowSums(raw.RC.filter1[,c(13:24)] >= 1, na.rm = T) >= 1
Rkeep <- rowSums(raw.RC.filter1[,c(25:36)] >= 1, na.rm = T) >= 1

# merge three logical vectors into one column by cbind(## & ##)
# cbind(## , ##) will combine several columns into one dataframe. 
FLRkeep <- cbind(Fkeep & Lkeep & Rkeep) 

raw.RC.filter2 <- raw.RC.filter1[FLRkeep,] 
dim(raw.RC.filter2)
message("Number of genes BEFORE filter2: ", dim(raw.RC.filter1)[1], 
        "\nNumber of genes AFTER filter2: ",  dim(raw.RC.filter2)[1])

# Create a new dataframe for the average value of each species/populations
head(raw.RC.filter2)
MeanRC <- data.frame("gene" = row.names(raw.RC.filter2))
head(MeanRC)

MeanRC$CG_F <- rowMeans(raw.RC.filter2[,CG_F], na.rm = T)
MeanRC$CR_F <- rowMeans(raw.RC.filter2[,CR_F], na.rm = T)
MeanRC$CO_F <- rowMeans(raw.RC.filter2[,CO_F], na.rm = T)
MeanRC$CG_L <- rowMeans(raw.RC.filter2[,CG_L], na.rm = T)
MeanRC$CR_L <- rowMeans(raw.RC.filter2[,CR_L], na.rm = T)
MeanRC$CO_L <- rowMeans(raw.RC.filter2[,CO_L], na.rm = T)
MeanRC$CG_R <- rowMeans(raw.RC.filter2[,CG_R], na.rm = T)
MeanRC$CR_R <- rowMeans(raw.RC.filter2[,CR_R], na.rm = T)
MeanRC$CO_R <- rowMeans(raw.RC.filter2[,CO_R], na.rm = T)

MeanRC$ASI_Co_F <- rowMeans(raw.RC.filter2[,ASI_Co_F], na.rm = T) 
MeanRC$EUR_Co_F <- rowMeans(raw.RC.filter2[,EUR_Co_F], na.rm = T) 
MeanRC$ME_Co_F <- rowMeans(raw.RC.filter2[,ME_Co_F], na.rm = T) 
MeanRC$CASI_Co_F <- rowMeans(raw.RC.filter2[,CASI_Co_F], na.rm = T) 
MeanRC$ASI_Cg_F <- rowMeans(raw.RC.filter2[,ASI_Cg_F], na.rm = T) 
MeanRC$EUR_Cg_F <- rowMeans(raw.RC.filter2[,EUR_Cg_F], na.rm = T) 
MeanRC$ME_Cg_F <- rowMeans(raw.RC.filter2[,ME_Cg_F], na.rm = T) 
MeanRC$CASI_Cg_F <- rowMeans(raw.RC.filter2[,CASI_Cg_F], na.rm = T) 

MeanRC$ASI_Co_L <- rowMeans(raw.RC.filter2[,ASI_Co_L], na.rm = T) 
MeanRC$EUR_Co_L <- rowMeans(raw.RC.filter2[,EUR_Co_L], na.rm = T) 
MeanRC$ME_Co_L <- rowMeans(raw.RC.filter2[,ME_Co_L], na.rm = T) 
MeanRC$CASI_Co_L <- rowMeans(raw.RC.filter2[,CASI_Co_L], na.rm = T) 
MeanRC$ASI_Cg_L <- rowMeans(raw.RC.filter2[,ASI_Cg_L], na.rm = T) 
MeanRC$EUR_Cg_L <- rowMeans(raw.RC.filter2[,EUR_Cg_L], na.rm = T) 
MeanRC$ME_Cg_L <- rowMeans(raw.RC.filter2[,ME_Cg_L], na.rm = T) 
MeanRC$CASI_Cg_L <- rowMeans(raw.RC.filter2[,CASI_Cg_L], na.rm = T) 

MeanRC$ASI_Co_R <- rowMeans(raw.RC.filter2[,ASI_Co_R], na.rm = T) 
MeanRC$EUR_Co_R <- rowMeans(raw.RC.filter2[,EUR_Co_R], na.rm = T) 
MeanRC$ME_Co_R <- rowMeans(raw.RC.filter2[,ME_Co_R], na.rm = T) 
MeanRC$CASI_Co_R <- rowMeans(raw.RC.filter2[,CASI_Co_R], na.rm = T) 
MeanRC$ASI_Cg_R <- rowMeans(raw.RC.filter2[,ASI_Cg_R], na.rm = T) 
MeanRC$EUR_Cg_R <- rowMeans(raw.RC.filter2[,EUR_Cg_R], na.rm = T) 
MeanRC$ME_Cg_R <- rowMeans(raw.RC.filter2[,ME_Cg_R], na.rm = T) 
MeanRC$CASI_Cg_R <- rowMeans(raw.RC.filter2[,CASI_Cg_R], na.rm = T) 
row.names(MeanRC) <- MeanRC$gene
MeanRC$gene <- NULL
head(MeanRC)
dim(MeanRC)
dim(na.omit(MeanRC))

Tissue <- c("F","L","R")
for (i in Tissue) {
  reads <- MeanRC[, which(sub(".*_","",colnames(MeanRC)) == i)]
  reads <- na.omit(reads)
  myGroup <- c("CG","CR","CO",rep("Cbp",8))
  myGroup <- factor(myGroup, levels = c("CG","CR","CO","Cbp"))
  message("my group in ", i, " is: ", myGroup)
  myDE <- DGEList(counts = reads, group = myGroup)
  # Normalizing the data
  myDE <- calcNormFactors(myDE, method = "TMM")
  # Estimating the Dispersion
  # This is the first major step in the analysis of DGE data
  myDE <- estimateCommonDisp(myDE)
  # Extract TMM values.
  TMM <- myDE$pseudo.counts
  # Output data
  write.table(TMM, paste("InputData/AverageTMM_", i, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
  assign(paste0("TMM_", i), TMM)
}

# Combine TMM across all tissues.
TMM <- cbind(TMM_F, TMM_L, TMM_R)
head(TMM)
write.table(TMM, file = "InputData/Capsella_averageTMM_FLR.txt", quote = F, sep = "\t", row.names = T, col.names = T)

