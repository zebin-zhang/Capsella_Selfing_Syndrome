######################################################
#  pairwise comparison of CPM among Cg, Cr, and Co.  #
######################################################

## Set work directory and load files
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
#### This file contains reads counts infromation of each gene in all samples
raw.RC <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE_rename.txt", header=TRUE, row.names="gene")
head(raw.RC)
dim(raw.RC)
# CPMnorm <- read.table("Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv", header = T, sep = "\t")
# head(CPMnorm)
# dim(CPMnorm)
##### This file contains information of population and tissue in each sample
Accession <- read.csv("InputData/28genomes_annot_R.txt", header = T, sep = "\t")
head(Accession)

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
#####################################
## Filter the count data
# Filter1: Only retain genes that have at least one expression values for each populations

geneMis <- 1 # For a given gene, set the minimum samples number with expression values in each population / subpopulation
geneKeep <- cbind( (rowSums(raw.RC[,CG_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CR_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CO_F] >= 1, na.rm = T) >= geneMis) &
                   
                   (rowSums(raw.RC[,CG_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CR_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CO_L] >= 1, na.rm = T) >= geneMis) &
                   
                   (rowSums(raw.RC[,CG_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CR_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CO_R] >= 1, na.rm = T) >= geneMis) &
                     
                   (rowSums(raw.RC[,ASI_Co_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Co_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Co_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Co_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Cg_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Cg_F] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
                     
                   (rowSums(raw.RC[,ASI_Co_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Co_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Co_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Co_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Cg_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Cg_L] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
                     
                   (rowSums(raw.RC[,ASI_Co_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Co_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Co_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Co_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ASI_Cg_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,EUR_Cg_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,ME_Cg_R] >= 1, na.rm = T) >= geneMis) &
                   (rowSums(raw.RC[,CASI_Cg_R] >= 1, na.rm = T) >= geneMis) ) 
raw.RC.filter1 <- raw.RC[geneKeep,]
message(paste("Number of genes before filter1: ", dim(raw.RC)[1], "\nNumber of genes after filter1: ", dim(raw.RC.filter1)[1]))

head(raw.RC.filter1)
dim(raw.RC.filter1)

## Filter HSE
# Filter by NA number
# For each gene in each tissue in parental species, at least 5 samples contain none-zero value.
Fkeep <- rowSums(raw.RC.filter1[,c(1:16)] >= 1, na.rm = T) >= 5
Lkeep <- rowSums(raw.RC.filter1[,c(17:33)] >= 1, na.rm = T) >= 5
Rkeep <- rowSums(raw.RC.filter1[,c(34:48)] >= 1, na.rm = T) >= 5

# merge three logical vectors into one column by cbind(## & ##)
# cbind(## , ##) will combine several columns into one dataframe. 
FLRkeep <- cbind(Fkeep & Lkeep & Rkeep) 

raw.RC.filter2 <- raw.RC.filter1[FLRkeep,] 
dim(raw.RC.filter2)
message("Number of genes BEFORE filter2: ", dim(raw.RC.filter1), 
        "\nNumber of genes AFTER filter2: ",  dim(raw.RC.filter2))

### Calculate the library sizes 
# Library size for Cbp = (Co + Cg)*0.5:
libColsums <- colSums(raw.RC.filter2 + 1, na.rm = T) # add one because later we will have to log values.
data.lib.size <- c()
evenN <- seq(1,132, 2)
for (i in c(1:c(length(libColsums)))){
  if (i <= 36){
    data.lib.size[i] <- libColsums[i]
  }else{
    if (i %% 2 == 1) {
      CbpLib <- libColsums[i]+libColsums[i+1]
      data.lib.size[i] <- CbpLib*0.5
      data.lib.size[i+1] <- CbpLib*0.5
    }
  }
}
data.lib.size[length(libColsums)] <- (libColsums[i]+libColsums[i-1])*0.5 # for the last column
barplot(data.lib.size/1000000, names.arg=Accession$NewName, las =2, border = F,
        col= as.character(Accession$col), ylab = "Library size (million reads)", cex.names=0.5)
# pdf("Lib-size_norm_FRL_All2.pdf", height=2, width=5)
# par(mar=c(5, 4, 1, 0), cex=0.5)
# barplot(data.lib.size/1000000, names.arg=Accession$accession, las =2, border = F,
#      col= tissuesCol, ylab = "Library size (million reads)", cex.names=0.7)
# dev.off()


### Count per million (CPM) normilazation
CPMnorm <- raw.RC.filter2
logCPMnorm <- raw.RC.filter2
for (i in c(1:c(length(data.lib.size)))){
  message(data.lib.size[i])
  CPMnorm[,i] <- raw.RC.filter2[,i] / (data.lib.size[i] / 1000000) # why needs to +1 in here?
  logCPMnorm[,i] <- log((raw.RC.filter2[,i] / (data.lib.size[i] / 1000000)) + 1) # why needs to +1 in here?
}
head(CPMnorm)
dim(CPMnorm)
head(logCPMnorm)
dim(logCPMnorm)


write.table(CPMnorm, "InputData/Capsella_ReadsCount_CPM_Filtered.txt",  quote = F, sep = "\t", row.names = T)
write.table(logCPMnorm, "InputData/Capsella_ReadsCount_logCPM_Filtered.txt",  quote = F, sep = "\t", row.names = T)

popLablesHclust <- c("CR","CR","CR","CR","CG","CG","CG","CG","CO","CO","CO","CO",
                     "CR","CR","CR","CR","CG","CG","CG","CG","CO","CO","CO","CO",
                     "CR","CR","CR","CR","CG","CG","CG","CG","CO","CO","CO","CO",
                     "ME","ME","ASI","ASI","CASI","CASI","EUR","EUR","ME","ME","ME",
                     "ME","ASI","ASI","CASI","CASI","CASI","CASI","ASI","ASI","EUR",
                     "EUR","EUR","EUR","EUR","EUR","CASI","CASI","ME","ME","ASI",
                     "ASI","ME","ME","ASI","ASI","CASI","CASI","EUR","EUR","ME","ME",
                     "ME","ME","ASI","ASI","CASI","CASI","CASI","CASI","ASI","ASI",
                     "EUR","EUR","EUR","EUR","EUR","EUR","CASI","CASI","ME","ME",
                     "ASI","ASI","ME","ME","ASI","ASI","CASI","CASI","EUR","EUR",
                     "ME","ME","ME","ME","ASI","ASI","CASI","CASI","CASI","CASI",
                     "ASI","ASI","EUR","EUR","EUR","EUR","EUR","EUR","CASI","CASI","ME","ME","ASI","ASI")
specLabelsHclust <- as.character(Accession$species)

paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep="")

dim(CPMnorm)

# install.packages("miscTools")
library(miscTools)
colMedians(CPMnorm, na.rm = T)

barplot(colMeans(CPMnorm, na.rm = T),  
        las =2, 
        border = F,
        col= as.character(Accession$col), ylab = "mean CPM", cex.names=0.5,
        main = "Mean Read counts CPM normalise",
        names.arg=paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep=""))

barplot(colMedians(CPMnorm, na.rm = T),  
        las =2, 
        border = F,
        col= as.character(Accession$col), ylab = "Median CPM", cex.names=0.5,
        main = "Median Read counts CPM normalise",
        names.arg=paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep=""))

barplot(log(colMeans(CPMnorm, na.rm = T)),  las =2, border = F,
        col= as.character(Accession$col), ylab = "mean logCPM", cex.names=0.5,
        main = "Mean Read counts CPM normalise log transformation",
        names.arg=paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep=""))


barplot(log(colMedians(CPMnorm, na.rm = T)),  las =2, border = F,
        col= as.character(Accession$col), ylab = "Median logCPM", cex.names=0.5,
        main = "Median Read counts CPM normalise log transformation",
        names.arg=paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep=""))
##############################
# Define populations
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
##############################

CPMnorm.Pop <- data.frame(Gene=row.names(CPMnorm))
head(CPMnorm.Pop)
# Use the MEAN value of each sample to represents as population value
# Flower of Cr, Cg, Co
CPMnorm.Pop$CR_F <- rowMeans(CPMnorm[,CR_F], na.rm = T)
CPMnorm.Pop$CG_F <- rowMeans(CPMnorm[,CG_F], na.rm = T)
CPMnorm.Pop$CO_F <- rowMeans(CPMnorm[,CO_F], na.rm = T)
# Leaf of Cr, Cg, Co
CPMnorm.Pop$CR_L <- rowMeans(CPMnorm[,CR_L], na.rm = T)
CPMnorm.Pop$CG_L <- rowMeans(CPMnorm[,CG_L], na.rm = T)
CPMnorm.Pop$CO_L <- rowMeans(CPMnorm[,CO_L], na.rm = T)
# Root of Cr, Cg, Co
CPMnorm.Pop$CR_R <- rowMeans(CPMnorm[,CR_R], na.rm = T)
CPMnorm.Pop$CG_R <- rowMeans(CPMnorm[,CG_R], na.rm = T)
CPMnorm.Pop$CO_R <- rowMeans(CPMnorm[,CO_R], na.rm = T)
# Flower of Cbp_co
CPMnorm.Pop$Cbp_ASI_Co_F <- rowMeans(CPMnorm[,ASI_Co_F], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Co_F <- rowMeans(CPMnorm[,EUR_Co_F], na.rm = T)
CPMnorm.Pop$Cbp_ME_Co_F <- rowMeans(CPMnorm[,ME_Co_F], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Co_F <- rowMeans(CPMnorm[,CASI_Co_F], na.rm = T)
# Flower of Cbp_cg
CPMnorm.Pop$Cbp_ASI_Cg_F <- rowMeans(CPMnorm[,ASI_Cg_F], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Cg_F <- rowMeans(CPMnorm[,EUR_Cg_F], na.rm = T)
CPMnorm.Pop$Cbp_ME_Cg_F <- rowMeans(CPMnorm[,ME_Cg_F], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Cg_F <- rowMeans(CPMnorm[,CASI_Cg_F], na.rm = T)
# Leaf of Cbp_co
CPMnorm.Pop$Cbp_ASI_Co_L <- rowMeans(CPMnorm[,ASI_Co_L], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Co_L <- rowMeans(CPMnorm[,EUR_Co_L], na.rm = T)
CPMnorm.Pop$Cbp_ME_Co_L <- rowMeans(CPMnorm[,ME_Co_L], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Co_L <- rowMeans(CPMnorm[,CASI_Co_L], na.rm = T)
# Leaf of Cbp_cg
CPMnorm.Pop$Cbp_ASI_Cg_L <- rowMeans(CPMnorm[,ASI_Cg_L], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Cg_L <- rowMeans(CPMnorm[,EUR_Cg_L], na.rm = T)
CPMnorm.Pop$Cbp_ME_Cg_L <- rowMeans(CPMnorm[,ME_Cg_L], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Cg_L <- rowMeans(CPMnorm[,CASI_Cg_L], na.rm = T)
# Root of Cbp_co
CPMnorm.Pop$Cbp_ASI_Co_R <- rowMeans(CPMnorm[,ASI_Co_R], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Co_R <- rowMeans(CPMnorm[,EUR_Co_R], na.rm = T)
CPMnorm.Pop$Cbp_ME_Co_R <- rowMeans(CPMnorm[,ME_Co_R], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Co_R <- rowMeans(CPMnorm[,CASI_Co_R], na.rm = T)
# Root of Cbp_cg
CPMnorm.Pop$Cbp_ASI_Cg_R <- rowMeans(CPMnorm[,ASI_Cg_R], na.rm = T)
CPMnorm.Pop$Cbp_EUR_Cg_R <- rowMeans(CPMnorm[,EUR_Cg_R], na.rm = T)
CPMnorm.Pop$Cbp_ME_Cg_R <- rowMeans(CPMnorm[,ME_Cg_R], na.rm = T)
CPMnorm.Pop$Cbp_CASI_Cg_R <- rowMeans(CPMnorm[,CASI_Cg_R], na.rm = T)
row.names(CPMnorm.Pop) <- as.character(CPMnorm.Pop$Gene)
CPMnorm.Pop$Gene <- NULL
head(CPMnorm.Pop)
dim(CPMnorm.Pop)
dim(na.omit(CPMnorm.Pop))
str(CPMnorm.Pop)

## Comparasion 
CPMnorm.Pop.Comp <- CPMnorm.Pop[,c("CR_F","CG_F","CO_F","CR_L","CG_L","CO_L","CR_R","CG_R","CO_R")]
# First calcuate the absolute value cross species.
### Flower
CPMnorm.Pop.Comp$abs.F_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_F"] - CPMnorm.Pop.Comp[,"CO_F"])
CPMnorm.Pop.Comp$abs.F_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_F"] - CPMnorm.Pop.Comp[,"CR_F"])
CPMnorm.Pop.Comp$abs.F_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_F"] - CPMnorm.Pop.Comp[,"CO_F"])
### Leaf
CPMnorm.Pop.Comp$abs.L_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_L"] - CPMnorm.Pop.Comp[,"CO_L"])
CPMnorm.Pop.Comp$abs.L_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_L"] - CPMnorm.Pop.Comp[,"CR_L"])
CPMnorm.Pop.Comp$abs.L_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_L"] - CPMnorm.Pop.Comp[,"CO_L"])
### Root
CPMnorm.Pop.Comp$abs.R_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_R"] - CPMnorm.Pop.Comp[,"CO_R"])
CPMnorm.Pop.Comp$abs.R_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_R"] - CPMnorm.Pop.Comp[,"CR_R"])
CPMnorm.Pop.Comp$abs.R_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_R"] - CPMnorm.Pop.Comp[,"CO_R"])

## Pick the samllar value of CgCr & CgCo, as the differences from outcorsser to selfer
## Use the value of CrCo as the difference within selfer
# install.packages("matrixStats")
library(matrixStats)
# Flower
CPMnorm.Pop.Comp$F.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.F_CgCo","abs.F_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$F.sf2sf <- CPMnorm.Pop.Comp$abs.F_CrCo
# Leaf
CPMnorm.Pop.Comp$L.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.L_CgCo","abs.L_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$L.sf2sf <- CPMnorm.Pop.Comp$abs.L_CrCo
# Root
CPMnorm.Pop.Comp$R.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.R_CgCo","abs.R_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$R.sf2sf <- CPMnorm.Pop.Comp$abs.R_CrCo
# Since we want to find some genes simliar expressed within selfer but different from selfer to outcrosser.
# So we require the value of sf2sf the smaller the better and value of oc2sf the bigger the better
# Here I creat a index of differences between two mating system:
# D.mate =  oc2sf - sf2sf
CPMnorm.Pop.Comp$F_D.mate <- CPMnorm.Pop.Comp$F.oc2sf - CPMnorm.Pop.Comp$F.sf2sf
CPMnorm.Pop.Comp$L_D.mate <- CPMnorm.Pop.Comp$L.oc2sf - CPMnorm.Pop.Comp$L.sf2sf
CPMnorm.Pop.Comp$R_D.mate <- CPMnorm.Pop.Comp$R.oc2sf - CPMnorm.Pop.Comp$R.sf2sf
# Considering the reads count mapped to each gene are different, it's hard to speak the same D.mate will have the same effect on
# mating type transformation from outcrosser to selfer. For instance:
# gene1  Cg = 10  Cr = 5  Co = 6, so the oc2sf=10-6=4, sf2sf=6-5=1, the D.mate = 4-1 =3
# gene2  Cg = 1000  Cr = 995  Co = 996, so the oc2sf=1000-996=4, sf2sf=996-995=1, the D.mate = 4-1 =3
# For gene1, it does have differences between outcrosser and selfer, but not gene2.
# So, to minimize the bias form reads count, we shink the D value to from -1 to 1 by devide to the Max read count of that gene across species.
# D = D.mate / CPM.species,  where CPM.sepcies = Max(CPM.Cg, CPM.Cr, CPM.Co)
# When D = 0, means totally nothing difference between outcrosser vs. selfer and selfer vs. selfer;
# when 0 < D <= 1, means the difference of outcrosser vs. selfer is bigger than the difference of selfer vs. selfer;
# when -1 <= D < 0, means the difference within selfer is bigger than differences between selfer and outcrosser. 

# Maximum CPM amongst thress species per tissue.
CPMnorm.Pop.Comp$F_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CR_F","CG_F","CO_F")]), na.rm = TRUE)
CPMnorm.Pop.Comp$L_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CR_L","CG_L","CO_L")]), na.rm = TRUE)
CPMnorm.Pop.Comp$R_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CR_R","CG_R","CO_R")]), na.rm = TRUE)

# Perform the calculation
### Flower
CPMnorm.Pop.Comp$F.D <- CPMnorm.Pop.Comp$F_D.mate / CPMnorm.Pop.Comp$F_CPM.species
CPMnorm.Pop.Comp$L.D <- CPMnorm.Pop.Comp$L_D.mate / CPMnorm.Pop.Comp$L_CPM.species
CPMnorm.Pop.Comp$R.D <- CPMnorm.Pop.Comp$R_D.mate / CPMnorm.Pop.Comp$R_CPM.species

# HIST
hist(CPMnorm.Pop.Comp$F.D, breaks = 100)
hist(CPMnorm.Pop.Comp$L.D, breaks = 100)
hist(CPMnorm.Pop.Comp$R.D, breaks = 100)

head(CPMnorm.Pop.Comp)
# Statistic how many genes with D value bigger than 0
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,])[1]
message("In Flower, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,])[1])
message("In Leaf, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,])[1])
message("In Root, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,])[1])

# Pick the top 5% gene
library(ggplot2)
library(magrittr)
library(dplyr)
#library(ps)
#install.packages("hrbrthemes")
#library(hrbrthemes)
CPMnorm.Pop.Comp.Flower.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(F.D > 0) %>%
  subset(F.D > quantile(F.D, prob = 1 - 5/100))
CPMnorm.Pop.Comp.Leaf.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(L.D > 0) %>%
  subset(L.D > quantile(L.D, prob = 1 - 5/100))
CPMnorm.Pop.Comp.Root.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(R.D > 0) %>%
  subset(R.D > quantile(R.D, prob = 1 - 5/100))
# Pick the cutoff value
min(CPMnorm.Pop.Comp.Flower.top0.05$F.D)  
min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D)  
min(CPMnorm.Pop.Comp.Root.top0.05$R.D)  
message("In flower, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Flower.top0.05$F.D))
message("In leaf, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D))
message("In root, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Root.top0.05$R.D))

# Group genes by D index
CPMnorm.Pop.Comp$F.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,]$F.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$F.D > min(CPMnorm.Pop.Comp.Flower.top0.05$F.D)),]$F.D.group <- "Top 5%"
CPMnorm.Pop.Comp$F.D.group <- factor(CPMnorm.Pop.Comp$F.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$F.D.group)

CPMnorm.Pop.Comp$L.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,]$L.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$L.D > min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D)),]$L.D.group <- "Top 5%"
CPMnorm.Pop.Comp$L.D.group <- factor(CPMnorm.Pop.Comp$L.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$L.D.group)

CPMnorm.Pop.Comp$R.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,]$R.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$R.D > min(CPMnorm.Pop.Comp.Root.top0.05$R.D)),]$R.D.group <- "Top 5%"
CPMnorm.Pop.Comp$R.D.group <- factor(CPMnorm.Pop.Comp$R.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$R.D.group)

head(CPMnorm.Pop.Comp)
mean(CPMnorm.Pop.Comp$F.D)
median(CPMnorm.Pop.Comp$F.D)
mean(CPMnorm.Pop.Comp$L.D)
median(CPMnorm.Pop.Comp$L.D)
mean(CPMnorm.Pop.Comp$R.D)
median(CPMnorm.Pop.Comp$R.D)


p1 <- ggplot(CPMnorm.Pop.Comp, aes(x=F.D, color = "black", fill = F.D.group)) +
   geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
   scale_fill_manual(values=c("red", "lightpink","gray90")) +
   labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
        y = "Gene counts", title = "Flower") +
   labs(fill="") + # remove the title for legend
   theme_bw() +
   theme(plot.title = element_text(hjust = 0.5), # put title to centre
         legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
p2 <- ggplot(CPMnorm.Pop.Comp, aes(x=L.D, color = "black", fill = L.D.group)) +
   geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
   scale_fill_manual(values=c("green4", "#D7FFD7","gray90")) +
   labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
        y = "Gene counts", title = "Leaf") +
   labs(fill="") + # remove the title for legend
   theme_bw() +
   theme(plot.title = element_text(hjust = 0.5), # put title to centre
         legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
p3 <- ggplot(CPMnorm.Pop.Comp, aes(x=R.D, color = "black", fill = R.D.group)) +
   geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
   scale_fill_manual(values=c("purple4", "#EACDFE","gray90")) +
   labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
        y = "Gene counts", title = "Root") +
   labs(fill="") + # remove the title for legend
   theme_bw() +
   theme(plot.title = element_text(hjust = 0.5), # put title to centre
         legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
library("grid")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p3, vp = vplayout(1,3))

# Extract gene names in top 5% regions
# In flower
head(CPMnorm.Pop.Comp.Flower.top0.05)
Top.5.genes.Flower <- row.names(CPMnorm.Pop.Comp.Flower.top0.05)
CPMnorm.Pop.Comp.Flower.top0.05$Gene.name <- rownames(CPMnorm.Pop.Comp.Flower.top0.05)
dim(CPMnorm.Pop.Comp.Flower.top0.05)
# In Leaf
head(CPMnorm.Pop.Comp.Leaf.top0.05)
Top.5.genes.Leaf <- row.names(CPMnorm.Pop.Comp.Leaf.top0.05)
CPMnorm.Pop.Comp.Leaf.top0.05$Gene.name <- rownames(CPMnorm.Pop.Comp.Leaf.top0.05)
dim(CPMnorm.Pop.Comp.Leaf.top0.05)
# In Root
head(CPMnorm.Pop.Comp.Root.top0.05)
Top.5.genes.Root <- row.names(CPMnorm.Pop.Comp.Root.top0.05)
CPMnorm.Pop.Comp.Root.top0.05$Gene.name <- rownames(CPMnorm.Pop.Comp.Root.top0.05)
dim(CPMnorm.Pop.Comp.Root.top0.05)

message("The gene number locate in the top 5% regions of flower are: ", dim(CPMnorm.Pop.Comp.Flower.top0.05)[1])
message("The gene number locate in the top 5% regions of leaf are: ", dim(CPMnorm.Pop.Comp.Leaf.top0.05)[1])
message("The gene number locate in the top 5% regions of root are: ", dim(CPMnorm.Pop.Comp.Root.top0.05)[1])
# Convert gene names of Capsella to Arabidopsis thaliana 拟南芥
cr2at <- read.table("AT_CR_mart_export.txt", header=TRUE, sep = "\t")

cr2at <- cr2at[,c("Ortholog.gene.name", "Ortholog.organism.name", "Gene.name", "Organism.name", "Relationship")]
names(cr2at) <- c("Gene.name", "Organism.name", "Ortholog.gene.name", "Ortholog.organism.name", "Relationship")
head(cr2at)
dim(cr2at)

Flower.geneKeep <- merge(CPMnorm.Pop.Comp.Flower.top0.05, cr2at, by = "Gene.name")
Leaf.geneKeep <- merge(CPMnorm.Pop.Comp.Leaf.top0.05, cr2at, by = "Gene.name")
Root.geneKeep <- merge(CPMnorm.Pop.Comp.Root.top0.05, cr2at, by = "Gene.name")

Flower.geneKeep <- Flower.geneKeep[,c(1, 35:41, 2:34)]
Leaf.geneKeep <- Leaf.geneKeep[,c(1, 35:41, 2:34)]
Root.geneKeep <- Root.geneKeep[,c(1, 35:41, 2:34)]
head(Flower.geneKeep)
dim(Flower.geneKeep)


message("In Flower, The gene counts within top 5% region is: ", dim(CPMnorm.Pop.Comp.Flower.top0.05)[1], 
        "; Of those, ", dim(Flower.geneKeep)[1], " genes have ortholog names on A. thaliana.")
message("In Leaf, The gene counts within top 5% region is: ", dim(CPMnorm.Pop.Comp.Leaf.top0.05)[1], 
        "; Of those, ", dim(Leaf.geneKeep)[1], " genes have ortholog names on A. thaliana.")
message("In Root, The gene counts within top 5% region is: ", dim(CPMnorm.Pop.Comp.Root.top0.05)[1], 
        "; Of those, ", dim(Root.geneKeep)[1], " genes have ortholog names on A. thaliana.")

# write.table(Flower.geneKeep, file = "Top_5%_differences_from_outcrosser_to_selfer_in_Flower.txt", quote = F, sep = "\t", row.names = F)
# write.table(Leaf.geneKeep, file = "Top_5%_differences_from_outcrosser_to_selfer_in_Leaf.txt", quote = F, sep = "\t", row.names = F)
# write.table(Root.geneKeep, file = "Top_5%_differences_from_outcrosser_to_selfer_in_Root.txt", quote = F, sep = "\t", row.names = F)

head(CPMnorm.Pop.Comp)
CPMnorm.Pop.Comp$Gene.name <- row.names(CPMnorm.Pop.Comp)
# write.table(CPMnorm.Pop.Comp, file = "OutputData/CPMnorm.Pop.Comp.txt", quote = F, sep = "\t", row.names = F)
# Heatmap cutoff 
AllTop5 <- CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$F.D.group == "Top 5%") | (CPMnorm.Pop.Comp$L.D.group == "Top 5%") | (CPMnorm.Pop.Comp$R.D.group == "Top 5%"), ]
AllTop5$Gene.name <- row.names(AllTop5)
All.geneKeep <- merge(AllTop5, cr2at, by = "Gene.name")
All.geneKeep <- All.geneKeep[,c(1, 35:41, 2:34)]
head(All.geneKeep)
dim(All.geneKeep)
str(All.geneKeep)

library(pheatmap)
# CPMnorm.Pop.fc$gene <- row.names(CPMnorm.Pop.fc)
# head(CPMnorm.Pop.fc)
# dim(CPMnorm.Pop.fc)

All.heatmap <- All.geneKeep[,c("F.D","L.D","R.D")]
head(All.heatmap)
dim(All.heatmap)
All.heatmap <- All.heatmap[!is.infinite(All.heatmap$F.D),]
All.heatmap_subset <- as.matrix(All.heatmap)
pheatmap(All.heatmap_subset)
my_tissue_col <- data.frame(Tissue = rep(c("Flower", "Leaf","Root")))
tissue_colors = list(
  Tissue = c(Flower="firebrick1", Leaf="forestgreen", Root="purple4"))

#my_sample_col <- data.frame(sample = rep(c("Flower", "Leaf","Root"), c(3,3,3)))
row.names(my_tissue_col) <- colnames(All.heatmap_subset)
# Color pattern
paletteFunc <- colorRampPalette(c("navyblue",  "white", 'red'));
palette     <- paletteFunc(100);

pheatmap(All.heatmap_subset, annotation_col = my_tissue_col, cluster_cols =  F, show_rownames=F, col=palette,
         annotation_colors = tissue_colors)
# dim(fc.heatmap_subset)
# head(fc.heatmap_subset)

CompKeep <- All.heatmap[((All.heatmap$F.D > 0.3630918) & (All.heatmap$L.D < 0) & (All.heatmap$R.D < 0)),]
dim(CompKeep)
pheatmap(CompKeep, annotation_col = my_tissue_col, cluster_cols =  F, show_rownames=F, col=palette,
         annotation_colors = tissue_colors)


# Venn Diagram
install.packages("VennDiagram")
library(VennDiagram)
Flower.GeneName <- Flower.geneKeep$Gene.name
Leaf.GeneName <- Leaf.geneKeep$Gene.name
Root.GeneName <- Root.geneKeep$Gene.name

venn.diagram(
  x = list(Flower.GeneName, Leaf.GeneName, Root.GeneName),
  category.names = c("Flower top 5%" , "Leaf top 5% " , "Root top 5%"),
  filename = 'test_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  compression = "lzw",
  lwd = 1,
  col=c("firebrick1", 'forestgreen', 'purple4'),
  fill = c(alpha("firebrick1",0.3), alpha('forestgreen',0.3), alpha('purple4',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("firebrick1", 'forestgreen', 'purple4'),
  rotation = 1
)
# Flower comparison of Pascal and EdgeR
FlowerKeep.EdgeR
venn.diagram(
  x = list(Flower.GeneName, FlowerKeep.EdgeR),
  category.names = c("Flower top 5%" , "EdgeR top 5% " ),
  filename = 'Flower_EdgeR_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  col=c("firebrick1", 'navy'),
  fill = c(alpha("firebrick1",0.3), alpha('navy',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("firebrick1", 'navy')
)

###########################
#### The purpose of the comparasion is to finde some genes that colse between Cr and Co, but distance of both Cg vs Cr, and Cg vs Co.
#### Specificly are, for the two values of Cg vs Cr, and Cg vs Co, the bigger the better.
#### For Cr vs Co, the smaller the better.

#### The fold changes should be the minimum value of ( Cg vs Cr, and Cg vs Co ) minus by (Cr vs Co)

#### Flower
head(CPMnorm.Pop.Comp)
# CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F_CrCo == 0,]
# head(CPMnorm.Pop)
# CPMnorm["Carubv10002622m.g",]
CPMnorm.Pop.Comp$F_fc <- NA
CPMnorm.Pop.Comp$F_log2fc <- NA
for(i in c(1:dim(CPMnorm.Pop.Comp)[1])){
  if (CPMnorm.Pop.Comp[i,]$rel.F_CgCo >= CPMnorm.Pop.Comp[i,]$rel.F_CgCr) {
    CPMnorm.Pop.Comp[i,]$F_fc <- CPMnorm.Pop.Comp[i,]$rel.F_CgCr / CPMnorm.Pop.Comp[i,]$rel.F_CrCo
    CPMnorm.Pop.Comp[i,]$F_log2fc <- log2(CPMnorm.Pop.Comp[i,]$rel.F_CgCr / CPMnorm.Pop.Comp[i,]$rel.F_CrCo)
  }else{
    CPMnorm.Pop.Comp[i,]$F_fc <- CPMnorm.Pop.Comp[i,]$rel.F_CgCo / CPMnorm.Pop.Comp[i,]$rel.F_CrCo
    CPMnorm.Pop.Comp[i,]$F_log2fc <- log2(CPMnorm.Pop.Comp[i,]$rel.F_CgCo / CPMnorm.Pop.Comp[i,]$rel.F_CrCo)
  }
}

hist(CPMnorm.Pop.Comp$F_fc, breaks = 100)
median(CPMnorm.Pop.Comp$F_fc)
hist(CPMnorm.Pop.Comp$F_log2fc, breaks = 100)
#### Leaf
CPMnorm.Pop.Comp$L_fc <- NA
CPMnorm.Pop.Comp$L_log2fc <- NA
for(i in c(1:dim(CPMnorm.Pop.Comp)[1])){
  if (CPMnorm.Pop.Comp[i,]$L_CgCo >= CPMnorm.Pop.Comp[i,]$L_CgCr) {
    CPMnorm.Pop.Comp[i,]$L_fc <- CPMnorm.Pop.Comp[i,]$L_CgCr / CPMnorm.Pop.Comp[i,]$L_CrCo
    CPMnorm.Pop.Comp[i,]$L_log2fc <- log2(CPMnorm.Pop.Comp[i,]$L_CgCr / CPMnorm.Pop.Comp[i,]$L_CrCo)
  }else{
    CPMnorm.Pop.Comp[i,]$L_fc <- CPMnorm.Pop.Comp[i,]$L_CgCo / CPMnorm.Pop.Comp[i,]$L_CrCo
    CPMnorm.Pop.Comp[i,]$L_log2fc <- log2(CPMnorm.Pop.Comp[i,]$L_CgCo / CPMnorm.Pop.Comp[i,]$L_CrCo)
  }
}

hist(CPMnorm.Pop.Comp$L_fc, breaks = 100)
hist(CPMnorm.Pop.Comp$L_log2fc, breaks = 100)


#### Root
CPMnorm.Pop.Comp$R_fc <- NA
CPMnorm.Pop.Comp$R_log2fc <- NA
for(i in c(1:dim(CPMnorm.Pop.Comp)[1])){
  if (CPMnorm.Pop.Comp[i,]$R_CgCo >= CPMnorm.Pop.Comp[i,]$R_CgCr) {
    CPMnorm.Pop.Comp[i,]$R_fc <- CPMnorm.Pop.Comp[i,]$R_CgCr / CPMnorm.Pop.Comp[i,]$R_CrCo
    CPMnorm.Pop.Comp[i,]$R_log2fc <- log2(CPMnorm.Pop.Comp[i,]$R_CgCr / CPMnorm.Pop.Comp[i,]$R_CrCo)
  }else{
    CPMnorm.Pop.Comp[i,]$R_fc <- CPMnorm.Pop.Comp[i,]$R_CgCo / CPMnorm.Pop.Comp[i,]$R_CrCo
    CPMnorm.Pop.Comp[i,]$R_log2fc <- log2(CPMnorm.Pop.Comp[i,]$R_CgCo / CPMnorm.Pop.Comp[i,]$R_CrCo)
  }
}

hist(CPMnorm.Pop.Comp$R_fc, breaks = 100)
hist(CPMnorm.Pop.Comp$R_log2fc, breaks = 100)

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}
col[which(CPMnorm.Pop.Comp$F_log2fc >= 2.322)] <- "firebrick1"

par(mfrow=c(3,2))
hist(CPMnorm.Pop.Comp$F_fc, breaks = 100, main="Fold changes of distances from outcrosser to selfer in Flower", xlab="CPMnorm")
hist(CPMnorm.Pop.Comp$F_log2fc, breaks = 100, 
     main="Fold changes of distances from outcrosser to selfer in Flower", xlab="log2(CPMnorm)")
hist(CPMnorm.Pop.Comp$L_fc, breaks = 100, main="Fold changes of distances from outcrosser to selfer in Leaf", xlab="CPMnorm")
hist(CPMnorm.Pop.Comp$L_log2fc, breaks = 100, main="Fold changes of distances from outcrosser to selfer in Leaf", xlab="log2(CPMnorm)")
hist(CPMnorm.Pop.Comp$R_fc, breaks = 100, main="Fold changes of distances from outcrosser to selfer in Root", xlab="CPMnorm")
hist(CPMnorm.Pop.Comp$R_log2fc, breaks = 100, main="Fold changes of distances from outcrosser to selfer in Root", xlab="log2(CPMnorm)")

par(resetPar())
par("mfrow") 
dev.off()

# fold change cut off : FC >= 5, log2FC >= 2.322
CPMnorm.Pop.fc <- CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$F_fc >= 5) | (CPMnorm.Pop.Comp$L_fc >= 5 ) | (CPMnorm.Pop.Comp$R_fc >= 5), ]

###################
# Load libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
library(DESeq2)
library(reshape)
library(ggplot2)
library(ggrepel)
# BiocManager::install("DEGreport")
library(DEGreport)
library(RColorBrewer)
# install.packages("pheatmap")
library(pheatmap)
CPMnorm.Pop.fc$gene <- row.names(CPMnorm.Pop.fc)
head(CPMnorm.Pop.fc)
dim(CPMnorm.Pop.fc)

fc.heatmap <- CPMnorm.Pop.fc[,c("F_log2fc","L_log2fc","R_log2fc")]
head(fc.heatmap)
dim(fc.heatmap)
fc.heatmap <- fc.heatmap[!is.infinite(fc.heatmap$F_log2fc),]
fc.heatmap_subset <- as.matrix(fc.heatmap)
pheatmap(fc.heatmap_subset)
my_tissue_col <- data.frame(Tissue = rep(c("Flower", "Leaf","Root")))
tissue_colors = list(
  Tissue = c(Flower="firebrick1", Leaf="forestgreen", Root="peru"))

#my_sample_col <- data.frame(sample = rep(c("Flower", "Leaf","Root"), c(3,3,3)))
row.names(my_tissue_col) <- colnames(fc.heatmap_subset)
# Color pattern
paletteFunc <- colorRampPalette(c("navyblue",  "white", 'red'));
palette     <- paletteFunc(100);

pheatmap(fc.heatmap_subset, annotation_col = my_tissue_col, cluster_cols =  F, show_rownames=F, col=palette,
         annotation_colors = tissue_colors)
dim(fc.heatmap_subset)
head(fc.heatmap_subset)

CompKeep <- fc.heatmap[((fc.heatmap$F_log2fc > 3.322) & (fc.heatmap$L_log2fc < 0) & (fc.heatmap$R_log2fc < 0)),]
dim(CompKeep)
pheatmap(CompKeep, annotation_col = my_tissue_col, cluster_cols =  F, show_rownames=F, col=palette,
         annotation_colors = tissue_colors)

head(CompKeep)
dim(CompKeep)
hist(CompKeep$F_log2fc, breaks = 50, main="", xlab="log2(CPMnorm)")

at2cr <- read.table("AT_CR_mart_export.txt", header=TRUE, sep = "\t")
head(at2cr)
dim(at2cr)

Flower.gene.keep <- row.names(CompKeep)

Flower.gene.keep.2AT <- at2cr[at2cr$Ortholog.gene.name %in% Flower.gene.keep,]
head(Flower.gene.keep.2AT)
dim(Flower.gene.keep.2AT)

CompKeep$Ortholog.gene.name <- row.names(CompKeep)
head(CompKeep)

Flower.gene.keep.at2cr <- merge(CompKeep, Flower.gene.keep.2AT, by = "Ortholog.gene.name")
Flower.gene.keep.at2cr <- Flower.gene.keep.at2cr[,c("Ortholog.organism.name","Ortholog.gene.name","F_log2fc","L_log2fc","R_log2fc",
                                                    "Organism.name", "Gene.name", "Relationship")]
names(Flower.gene.keep.at2cr) <- c("Organism.name","Gene.name","F_log2fc","L_log2fc","R_log2fc",
                                   "Ortholog.organism.name","Ortholog.gene.name","Relationship")
head(Flower.gene.keep.at2cr)
dim(Flower.gene.keep.at2cr)

write.table(Flower.gene.keep.at2cr, file = "Flower.gene.keep.CR2AT.txt", quote = F, sep = "\t", row.names = F)


#### Split to Flower
SigComp.F <- melt(CPMnorm.Pop.fc[,c("F_CgCo","F_CgCr","F_CrCo","F_fc", "F_log2fc", "gene")], id.vars =c("gene", "F_fc", "F_log2fc"))
colnames(SigComp.F) <- c("gene","FoldChanges", "log2FoldChanges", "Comparison","RelativeDistance")
SigComp.F$tissue <- "Flower"
SigComp.F$Comparison <- sub("F_","",SigComp.F$Comparison)
head(SigComp.F)

#### Split to Leaf
SigComp.L <- melt(CPMnorm.Pop.fc[,c("L_CgCo","L_CgCr","L_CrCo","L_fc","gene")], id.vars =c("gene", "L_fc"))
colnames(SigComp.L) <- c("gene","FoldChanges","Comparison","RelativeDistance")
SigComp.L$tissue <- "Leaf"
SigComp.L$Comparison <- sub("l_","",SigComp.L$Comparison)
head(SigComp.L)
#### Split to Root
SigComp.R <- melt(CPMnorm.Pop.fc[,c("R_CgCo","R_CgCr","R_CrCo","R_fc","gene")], id.vars =c("gene", "R_fc"))
colnames(SigComp.R) <- c("gene","FoldChanges","Comparison","RelativeDistance")
SigComp.R$tissue <- "Root"
SigComp.R$Comparison <- sub("R_","",SigComp.R$Comparison)
head(SigComp.R)

#### combine files by tissues
SigComp <- rbind(SigComp.F,SigComp.L,SigComp.R)
head(SigComp)




