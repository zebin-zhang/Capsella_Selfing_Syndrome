#########################################################################
#          Capsella Project Analysis Steps -- Zebin Zhang               #
#########################################################################

##  Set work directly and load files.
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project")
# This file contains reads counts infromation of each gene in all samples
RNAreadsCount <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv", header=TRUE, row.names="gene")
head(RNAreadsCount)
# Differences between the two subgenomes HSE (homeologue specific expression)
hseD <- read.table("ASE_sign_prop_merged_masked.csv", header=TRUE, row.names="gene")
head(hseD)
dim(hseD)
# This file contails informations of population and tissue in each accession
Accession <- read.csv("InputData/28genomes_annot_R.csv", header = T)
head(Accession)

head(RNAreadsCount)
dim(RNAreadsCount)

### Define populations
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

geneMis <- 1 # For a given gene, set the minmum samples number with expression values in each population / subpopulation
geneKeep <- cbind( (rowSums(!is.na(RNAreadsCount[,ASI_Co_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Co_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Co_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Co_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ASI_Cg_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Cg_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Cg_F]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Cg_F]) >= 1) >= geneMis) &
                     
                     (rowSums(!is.na(RNAreadsCount[,ASI_Co_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Co_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Co_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Co_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ASI_Cg_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Cg_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Cg_L]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Cg_L]) >= 1) >= geneMis) &
                     
                     (rowSums(!is.na(RNAreadsCount[,ASI_Co_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Co_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Co_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Co_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ASI_Cg_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,EUR_Cg_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,ME_Cg_R]) >= 1) >= geneMis) &
                     (rowSums(!is.na(RNAreadsCount[,CASI_Cg_R]) >= 1) >= geneMis))
RNAreadsCount.filter1 <- RNAreadsCount[geneKeep,]
message(paste("Number of genes before filter1: ", dim(RNAreadsCount)[1], "\nNumber of genes after filter1: ", dim(RNAreadsCount.filter1)[1]))

head(RNAreadsCount.filter1)
dim(RNAreadsCount.filter1)

#####################################
## Filter2:  filter genes with significant HSE in Cbp
# Look at the distribution of significant HSE in Cbp
genesKeep <- rownames(RNAreadsCount.filter1)
hseDF <- hseD[rownames(hseD) %in% genesKeep, ]
hist(rowSums(is.na(hseDF[c(1:16)]) >= 1), main = "Cbp Flowers")
hist(rowSums(is.na(hseDF[c(17:33)]) >= 1), main = "Cbp Leaves")
hist(rowSums(is.na(hseDF[c(34:48)]) >= 1), main = "Cbp Roots")

## Filter HSE
# Filter by NA number
# For each gene in each tissue, at least 5 samples contain none-zero value.
HSEkeepF <- rowSums(!is.na(hseDF[,c(1:16)]) >= 1) >= 5
HSEkeepL <- rowSums(!is.na(hseDF[,c(17:33)]) >= 1) >= 5
HSEkeepR <- rowSums(!is.na(hseDF[,c(34:48)]) >= 1) >= 5

hseF <- hseDF[HSEkeepF,c(1:16)]
hseL <- hseDF[HSEkeepL,c(17:33)]
hseR <- hseDF[HSEkeepR,c(34:48)]
hseFLR <- hseDF[genesKeep,] # this is the filter1, not filter2, not the same with combine hseF hseL hseR

message(paste("Number of genes before filter2: ", dim(hseDF)[1]))
message(paste("Number of genes after filter2: ", dim(hseF)[1]))
message(paste("Number of genes after filter2: ", dim(hseL)[1]))
message(paste("Number of genes after filter2: ", dim(hseR)[1]))
message(paste("Number of genes after filter2: ", dim(hseFLR)[1]))

## HSE tissues mean
head(hseF)
hseNF <- colSums(!(is.na(hseF)))
hseNL <- colSums(!(is.na(hseL)))
hseNR <- colSums(!(is.na(hseR)))
hseN <- as.data.frame(cbind((c(hseNF,hseNL,hseNR)), c(rep("F",16), rep("L",16), rep('R',16))))
mean(hseNF)
mean(hseNL)
mean(hseNR)

### ANOVA
hseN$V1 <- as.numeric(hseN$V1)
hseN$V2 <- as.factor(hseN$V2)
summary(aov(hseN$V1~hseN$V2))

# Extract gene names after the filter1.
dim(hseFLR)
hseG <- rownames(hseFLR)

#####################################

### Calculate the library sizes 
# Library size for Cbp = (Co + Cg)*0.5:
libColsums <- colSums(RNAreadsCount.filter1 + 1, na.rm = T) # add one because later we will have to log values.
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
barplot(data.lib.size/1000000, names.arg=Accession$accession, las =2, border = F,
        col= as.character(Accession$col), ylab = "Library size (million reads)", cex.names=0.5)
# pdf("Lib-size_norm_FRL_All2.pdf", height=2, width=5)
# par(mar=c(5, 4, 1, 0), cex=0.5)
# barplot(data.lib.size/1000000, names.arg=annot$accession, las =2, border = F,
#        col= tissuesCol, ylab = "Library size (million reads)", cex.names=0.7)
# dev.off()

### Count per million (CPM) normilazation
CPMnorm <- RNAreadsCount.filter1
for (i in c(1:c(length(data.lib.size)))){
  CPMnorm[,i] <- log(((RNAreadsCount.filter1[,i]+1)/data.lib.size[i])*1000000) # why needs to +1 in here?
}
head(CPMnorm)
dim(CPMnorm)

# write.table(CPMnorm, "Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv",  quote = F, sep = "\t", row.names = T)
paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep="")

dim(CPMnorm)

barplot(colMeans(CPMnorm, na.rm = T),  las =2, border = F,
        col= as.character(Accession$col), ylab = "mean log CPM", cex.names=0.5,
        names.arg=paste(specLabelsHclust, "_", Accession$tissue, "_", popLablesHclust, sep=""))

### Hierarchical Cluster Analysis -- Clustering in samples

clustSamp <- dist(t(CPMnorm), diag = T, method = "euclidean")
hSamp <- hclust(clustSamp, method = "average")
plot(hSamp, main = "", xlab = "", cex=0.4)

### Define HSE color
dotB <- rgb(0,0,0,0.1)
dotR <- rgb(1,0,0,0.1)

allGenes <- rownames(CPMnorm)
geneCol <- c()
for (g in allGenes){
  if (g %in% hseG){
    geneCol <- c(geneCol, dotR)
  }else{
    geneCol <- c(geneCol, dotB)
  }
}

### Clustering in samples
clustSamp <- dist(t(CPMnorm), diag = T, method = "euclidean")
hSamp <- hclust(clustSamp, method = "average")

plot(hSamp, main = "", xlab = "", cex=0.4)

# Write image as separate file
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

# pdf("norm_pvclust_FRL_spec_mean_Samp.pdf", height=6, width=12)
# par(mar=c(0, 4, 1, 0), cex=1)
plot(hSamp, main = "", xlab = "", labels=paste(specLabelsHclust, "_", annot$accession, "_", popLablesHclust, sep=""), cex=0.5)
# dev.off()

## Calculate mean value
dim(CPMnorm)
head(CPMnorm)
str(CPMnorm)
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

###############################
# PCA in all tissues: Flower, Leaf, Root
# Perform PCA
pca.all <- prcomp(t(na.omit(CPMnorm.Pop)))
# NOTE: By default, prcomp() expects the samples to be rows and the genes to be columns
# So the data have to be transpose the matrix using the t() function
# If don't transpose the matrix, you will ultimately get a graph that shows how the genes are related to each other.

# prcomp() returns thress things:
    # 1) x -- x contains the principal components (PCs) for drawing a graph.
    ## plot PC1 and PC2
plot(pca.all$x[,1],pca.all$x[,2])
    # 2) sdev -- "standard deviation", to calculate how much variation in the original data each PC accounts for.
    ## how much variation in the original data PC1 accounts for.
pca.all.var <- pca.all$sdev^2
head(pca.all.var)
    ## Since the percentage of variation that each PC accounts for is way more interesting than the acutal value,
    ## we calculate the percentage..
pca.all.var.per <- round(pca.all.var/sum(pca.all.var)*100,1)
    ## Plotting the percentage is easy with barplot()
barplot(pca.all.var.per, main = "Contribuation of PCs", xlab = "Pricipal Component", ylab = "Percent Variation")
    # 3) rotation -- the loading scores rotation

library(ggplot2)
# Format the data the way ggplot2 likes it
pca.all.data <- data.frame(PC1=pca.all$x[,1],
                           PC2=pca.all$x[,2])
pca.all.data$ID <- row.names(pca.all.data)
pca.all.data$tissue <- c(rep("Flower",3),rep("Leaf",3),rep("Root",3),rep("Flower",8),rep("Leaf",8),rep("Root",8))
pca.all.data$population <- c(rep(c("CR","CG","CO"),3),
                             rep(c("ASI_Co","EUR_Co","ME_Co","CASI_Co",
                                   "ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg"),3))
pca.all.data$population <- factor(pca.all.data$population, levels = c("CG","CR","CO",
                                                                      "ASI_Cg","CASI_Cg","EUR_Cg","ME_Cg",
                                                                      "ASI_Co","CASI_Co","EUR_Co","ME_Co"))
pca.all.data$species <- c(rep(c("CR","CG","CO"),3),
                          rep(c(rep("Cbp_Co",4), rep("Cbp_Cg",4)),3))

head(pca.all.data)
str(pca.all.data)


# graph parttern 1
ggplot(data = pca.all.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.all.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.all.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Flowers, Leaves, and Roots")
# graph parttern 2
ggplot(pca.all.data, aes(PC1, PC2, col = population, fill = population))  +
  geom_point(size = 3, aes(shape = population)) +
  scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
  scale_color_manual(values = c("indianred1","darkviolet","blue",
                                "green4","orange","red","royalblue",
                                "green4","orange","red","royalblue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue",
                               "green4","orange","red","royalblue",
                               "green4","orange","red","royalblue")) +
  stat_ellipse(aes(x= PC1, y=PC2, group = tissue), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.all.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.all.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Flowers, Leaves, and Roots")

# Using the loading scores to determine which genes have the largest effect on where samples are plotted in the PCA plot.
# Loading scores for each PC. 
# Loading scores for PC1, since it accounts for 67.8% of the variation in the data.
all.loading.scores <- pca.all$rotation[,1]
# the absoulte number of magnitude on each gene along PC1, 
# genes with large negative values push samples to the left, and large positive value push samples to the right.
all.gene_scores <- abs(all.loading.scores)
# sort the magnitudes of the loading scores, from high to low
all.gene_scores_ranked <- sort(all.gene_scores, decreasing = TRUE)
# The names of the top10 genes with the largest loading score magnitudes.
all.top_10_genes <- names(all.gene_scores_ranked[1:10])
all.top_10_genes
# show the scores (both positive and negative)
pca.all$rotation[all.top_10_genes,1]

###############################
## PCA in Flowers
### Define flower
head(CPMnorm.Pop)
Flower <- c("CR_F","CG_F","CO_F",
            "Cbp_ASI_Co_F","Cbp_EUR_Co_F","Cbp_ME_Co_F","Cbp_CASI_Co_F",
            "Cbp_ASI_Cg_F","Cbp_EUR_Cg_F","Cbp_ME_Cg_F","Cbp_CASI_Cg_F")
CPMnorm.Flower <- CPMnorm.Pop[,Flower]
head(CPMnorm.Flower)
dim(CPMnorm.Flower)
dim(na.omit(CPMnorm.Flower))
### perform PCA in Flower
pca.flower <- prcomp(t(na.omit(CPMnorm.Flower)))
# 1) X
## plot PC1 and PC2
plot(pca.flower$x[,1],pca.flower$x[,2])
# 2) sdev 
pca.flower.var <- pca.flower$sdev^2
head(pca.flower.var)
## calculate the percentage..
pca.flower.var.per <- round(pca.flower.var/sum(pca.flower.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.flower.var.per, main = "Contribuation of PCs in Flowers", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Flowers
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.flower$x)[1]
pca.flower.data <- data.frame(ID = 1:dim(pca.flower$x)[1])
for (i in 1:dim(pca.flower$x)[1]){
  pcs <- paste0("PC", i)
  pca.flower.data[[pcs]] <- pca.flower$x[,i]
}

pca.flower.data$ID <- row.names(pca.flower$x)

pca.flower.data$ID <- row.names(pca.flower.data)
pca.flower.data$tissue <- "Flower"
pca.flower.data$population <- c(c("CR","CG","CO"),
                             c("ASI_Co","EUR_Co","ME_Co","CASI_Co",
                               "ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg"))
pca.flower.data$population <- factor(pca.flower.data$population, levels = c("CG","CR","CO",
                                                                      "ASI_Cg","CASI_Cg","EUR_Cg","ME_Cg",
                                                                      "ASI_Co","CASI_Co","EUR_Co","ME_Co"))
pca.flower.data$species <- c(c("CR","CG","CO"),
                             rep("Cbp_Co",4),
                             rep("Cbp_Cg",4))


# graph parttern 1
ggplot(data = pca.flower.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.flower.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.flower.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Flowers, Leaves, and Roots")
# graph parttern 2
# PC1 vs PC2
ggplot(pca.flower.data, aes(PC1, PC2, col = population, fill = population))  +
  geom_point(size = 3, aes(shape = population)) +
  scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
  scale_color_manual(values = c("indianred1","darkviolet","blue",
                                "green4","orange","red","royalblue",
                                "green4","orange","red","royalblue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue",
                               "green4","orange","red","royalblue",
                               "green4","orange","red","royalblue")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.flower.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.flower.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Flowers")
# Other PCs
dim(pca.flowerr.data)
for (i in 1:(dim(pca.flowerr.data)[1] - 1)){
  pci <- paste0("PC", i)
  j <- i +1
  pcj <- paste0("PC", j)
  
  plot <- ggplot(pca.flowerr.data, aes(pca.flowerr.data[[pci]], pca.flowerr.data[[pcj]], col = population, fill = population))  +
    geom_point(size = 3, aes(shape = population)) +
    scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
    scale_color_manual(values = c("indianred1","darkviolet","blue",
                                  "green4","orange","red","royalblue",
                                  "green4","orange","red","royalblue"))  +
    scale_fill_manual(values = c("indianred1","darkviolet","blue",
                                 "green4","orange","red","royalblue",
                                 "green4","orange","red","royalblue")) +
    #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
    xlab(paste("PC", i," - ", pca.flowerr.var.per[i], "%", sep = "")) +
    ylab(paste("PC", j, " - ", pca.flowerr.var.per[j], "%", sep = "")) +
    theme_bw() +
    ggtitle("Gene expression in Flowers")
  print(plot)
}


### Using the loading scores to determine which genes have the largest effect on where samples are plotted in the Flower PCA plot.
# - Loading scores for each PC. 
# Loading scores for PC1, since it contribute the biggest of the variation in the data.
flower.loading.scores <- pca.flower$rotation[,1]
# the absoulte number of magnitude on each gene along PC1, 
# genes with large negative values push samples to the left, and large positive value push samples to the right.
flower.gene_scores <- abs(flower.loading.scores)
# sort the magnitudes of the loading scores, from high to low
flower.gene_scores_ranked <- sort(flower.gene_scores, decreasing = TRUE)
# The names of the top10 genes with the largest loading score magnitudes.
flower.top_10_genes <- names(flower.gene_scores_ranked[1:10])
flower.top_10_genes
# show the scores (both positive and negative)
pca.flower$rotation[flower.top_10_genes,1]

###############################
## PCA in Leaves
### Define leaf
head(CPMnorm.Pop)
Leaf <- c("CR_L","CG_L","CO_L",
          "Cbp_ASI_Co_L","Cbp_EUR_Co_L","Cbp_ME_Co_L","Cbp_CASI_Co_L",
          "Cbp_ASI_Cg_L","Cbp_EUR_Cg_L","Cbp_ME_Cg_L","Cbp_CASI_Cg_L")
CPMnorm.Leaf <- CPMnorm.Pop[,Leaf]
head(CPMnorm.Leaf)
dim(CPMnorm.Leaf)
dim(na.omit(CPMnorm.Leaf))
### perform PCA in Leaf
pca.leaf <- prcomp(t(na.omit(CPMnorm.Leaf)))
summary(pca.leaf)
# 1) X
## plot PC1 and PC2
plot(pca.leaf$x[,1],pca.leaf$x[,2])
# 2) sdev 
pca.leaf.var <- pca.leaf$sdev^2
head(pca.leaf.var)
## calculate the percentage..
pca.leaf.var.per <- round(pca.leaf.var/sum(pca.leaf.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.leaf.var.per, main = "Contribuation of PCs in Leaves", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Leaves
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.leaf$x)[1]
pca.leaf.data <- data.frame(ID = 1:dim(pca.leaf$x)[1])
for (i in 1:dim(pca.leaf$x)[1]){
  pcs <- paste0("PC", i)
  pca.leaf.data[[pcs]] <- pca.leaf$x[,i]
}

pca.leaf.data$ID <- row.names(pca.leaf$x)

pca.leaf.data$tissue <- "Leaf"
pca.leaf.data$population <- c(c("CR","CG","CO"),
                              c("ASI_Co","EUR_Co","ME_Co","CASI_Co",
                                "ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg"))
pca.leaf.data$population <- factor(pca.leaf.data$population, levels = c("CG","CR","CO",
                                                                        "ASI_Cg","CASI_Cg","EUR_Cg","ME_Cg",
                                                                        "ASI_Co","CASI_Co","EUR_Co","ME_Co"))
pca.leaf.data$species <- c(c("CR","CG","CO"),
                           rep("Cbp_Co",4),
                           rep("Cbp_Cg",4))


# graph parttern 1
ggplot(data = pca.leaf.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.leaf.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.leaf.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Leaves, Leaves, and Leaves")
# graph parttern 2
# PC1 vs PC2
ggplot(pca.leaf.data, aes(PC1, PC2, col = population, fill = population))  +
  geom_point(size = 3, aes(shape = population)) +
  scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
  scale_color_manual(values = c("indianred1","darkviolet","blue",
                                "green4","orange","red","royalblue",
                                "green4","orange","red","royalblue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue",
                               "green4","orange","red","royalblue",
                               "green4","orange","red","royalblue")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.leaf.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.leaf.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Leaves")
# Other PCs
dim(pca.leaf.data)
for (i in 1:(dim(pca.leaf.data)[1] - 1)){
  pci <- paste0("PC", i)
  j <- i +1
  pcj <- paste0("PC", j)
  
  plot <- ggplot(pca.leaf.data, aes(pca.leaf.data[[pci]], pca.leaf.data[[pcj]], col = population, fill = population))  +
    geom_point(size = 3, aes(shape = population)) +
    scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
    scale_color_manual(values = c("indianred1","darkviolet","blue",
                                  "green4","orange","red","royalblue",
                                  "green4","orange","red","royalblue"))  +
    scale_fill_manual(values = c("indianred1","darkviolet","blue",
                                 "green4","orange","red","royalblue",
                                 "green4","orange","red","royalblue")) +
    #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
    xlab(paste("PC", i," - ", pca.leaf.var.per[i], "%", sep = "")) +
    ylab(paste("PC", j, " - ", pca.leaf.var.per[j], "%", sep = "")) +
    theme_bw() +
    ggtitle("Gene expression in Leaves")
  print(plot)
}

### Using the loading scores to determine which genes have the largest effect on where samples are plotted in the Leaf PCA plot.
# - Loading scores for each PC. 
# Loading scores for PC1, since it contribute the biggest of the variation in the data.
leaf.loading.scores <- pca.leaf$rotation[,1]
# the absoulte number of magnitude on each gene along PC1, 
# genes with large negative values push samples to the left, and large positive value push samples to the right.
leaf.gene_scores <- abs(leaf.loading.scores)
# sort the magnitudes of the loading scores, from high to low
leaf.gene_scores_ranked <- sort(leaf.gene_scores, decreasing = TRUE)
# The names of the top10 genes with the largest loading score magnitudes.
leaf.top_10_genes <- names(leaf.gene_scores_ranked[1:10])
leaf.top_10_genes
# show the scores (both positive and negative)
pca.leaf$rotation[leaf.top_10_genes,1]



###############################
## PCA in Roots
### Define root
head(CPMnorm.Pop)
Root <- c("CR_R","CG_R","CO_R",
          "Cbp_ASI_Co_R","Cbp_EUR_Co_R","Cbp_ME_Co_R","Cbp_CASI_Co_R",
          "Cbp_ASI_Cg_R","Cbp_EUR_Cg_R","Cbp_ME_Cg_R","Cbp_CASI_Cg_R")
CPMnorm.Root <- CPMnorm.Pop[,Root]
head(CPMnorm.Root)
dim(CPMnorm.Root)
dim(na.omit(CPMnorm.Root))
### perform PCA in Root
pca.root <- prcomp(t(na.omit(CPMnorm.Root)))
# 1) X
## plot PC1 and PC2
plot(pca.root$x[,1],pca.root$x[,2])
# 2) sdev 
pca.root.var <- pca.root$sdev^2
head(pca.root.var)
## calculate the percentage..
pca.root.var.per <- round(pca.root.var/sum(pca.root.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.root.var.per, main = "Contribuation of PCs in Roots", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Roots
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.root$x)[1]
pca.root.data <- data.frame(ID = 1:dim(pca.root$x)[1])
for (i in 1:dim(pca.root$x)[1]){
  pcs <- paste0("PC", i)
  pca.root.data[[pcs]] <- pca.root$x[,i]
}

pca.root.data$ID <- row.names(pca.root$x)
pca.root.data$tissue <- "Root"
pca.root.data$population <- c(c("CR","CG","CO"),
                              c("ASI_Co","EUR_Co","ME_Co","CASI_Co",
                                "ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg"))
pca.root.data$population <- factor(pca.root.data$population, levels = c("CG","CR","CO",
                                                                        "ASI_Cg","CASI_Cg","EUR_Cg","ME_Cg",
                                                                        "ASI_Co","CASI_Co","EUR_Co","ME_Co"))
pca.root.data$species <- c(c("CR","CG","CO"),
                           rep("Cbp_Co",4),
                           rep("Cbp_Cg",4))


# graph parttern 1
ggplot(data = pca.root.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Roots")
# graph parttern 2
ggplot(pca.root.data, aes(PC1, PC2, col = population, fill = population))  +
  geom_point(size = 3, aes(shape = population)) +
  scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
  scale_color_manual(values = c("indianred1","darkviolet","blue",
                                "green4","orange","red","royalblue",
                                "green4","orange","red","royalblue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue",
                               "green4","orange","red","royalblue",
                               "green4","orange","red","royalblue")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Roots")

# Other PCs
dim(pca.root.data)
for (i in 1:(dim(pca.root.data)[1] - 1)){
  pci <- paste0("PC", i)
  j <- i +1
  pcj <- paste0("PC", j)
  
  plot <- ggplot(pca.root.data, aes(pca.root.data[[pci]], pca.root.data[[pcj]], col = population, fill = population))  +
    geom_point(size = 3, aes(shape = population)) +
    scale_shape_manual(values=c(8,10,4,22,23,21,24,0,5,1,2))+
    scale_color_manual(values = c("indianred1","darkviolet","blue",
                                  "green4","orange","red","royalblue",
                                  "green4","orange","red","royalblue"))  +
    scale_fill_manual(values = c("indianred1","darkviolet","blue",
                                 "green4","orange","red","royalblue",
                                 "green4","orange","red","royalblue")) +
    #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
    xlab(paste("PC", i," - ", pca.root.var.per[i], "%", sep = "")) +
    ylab(paste("PC", j, " - ", pca.root.var.per[j], "%", sep = "")) +
    theme_bw() +
    ggtitle("Gene expression in Roots")
  print(plot)
}


### Using the loading scores to determine which genes have the largest effect on where samples are plotted in the Root PCA plot.
# - Loading scores for each PC. 
# Loading scores for PC1, since it contribute the biggest of the variation in the data.
root.loading.scores <- pca.root$rotation[,1]
# the absoulte number of magnitude on each gene along PC1, 
# genes with large negative values push samples to the left, and large positive value push samples to the right.
root.gene_scores <- abs(root.loading.scores)
# sort the magnitudes of the loading scores, from high to low
root.gene_scores_ranked <- sort(root.gene_scores, decreasing = TRUE)
# The names of the top10 genes with the largest loading score magnitudes.
root.top_10_genes <- names(root.gene_scores_ranked[1:10])
root.top_10_genes
# show the scores (both positive and negative)
pca.root$rotation[root.top_10_genes,1]

