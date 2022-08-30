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
                    (rowSums(raw.RC[,CO_R] >= 1, na.rm = T) >= geneMis) 
)

raw.RC.filter1 <- raw.RC[geneKeep,]
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.filter1)[1]))

## Filter HSE
# Filter by NA number
# For each gene in each tissue in parental species, at least 5 samples contain none-zero value.
Fkeep <- rowSums(raw.RC.filter1[,c(1:12)] >= 1, na.rm = T) >= 5
Lkeep <- rowSums(raw.RC.filter1[,c(13:24)] >= 1, na.rm = T) >= 5
Rkeep <- rowSums(raw.RC.filter1[,c(25:36)] >= 1, na.rm = T) >= 5

# merge three logical vectors into one column by cbind(## & ##)
# cbind(## , ##) will combine several columns into one dataframe. 
FLRkeep <- cbind(Fkeep & Lkeep & Rkeep) 

raw.RC.filter2 <- raw.RC.filter1[FLRkeep,] 
dim(raw.RC.filter2)
message("Number of genes BEFORE filter2: ", dim(raw.RC.filter1)[1], 
        "\nNumber of genes AFTER filter2: ",  dim(raw.RC.filter2)[1])

ReadsKeep <- raw.RC.filter2
head(raw.RC.filter2)
DiplidKeep <- raw.RC.filter2[,1:36]
head(DiplidKeep)


Tissue <- c("F","L","R")
for (i in Tissue) {
  reads <- DiplidKeep[, which(sub(".*_","",colnames(DiplidKeep)) == i)]
  reads <- na.omit(reads)
  myGroup <- c(rep("CR",4), rep("CG",4), rep("CO",4) )
  myGroup <- factor(myGroup, levels = c("CG","CR","CO"))
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
  write.table(TMM, paste("InputData/Diploids_individual_TMM_", i, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
  assign(paste0("TMM_", i), TMM)
}

dim(TMM_F)
dim(TMM_L)
dim(TMM_R)

head(TMM_F)
TMM <- cbind(TMM_F, TMM_L, TMM_R)
head(TMM)

write.table(TMM, file = "InputData/Diploids_individual_TMM_FLR.txt", quote = F, sep = "\t", row.names = T, col.names = T)











