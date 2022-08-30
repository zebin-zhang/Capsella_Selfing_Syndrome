#########################################################################
#          Capsella Project Analysis Steps -- Zebin Zhang               #
#########################################################################

##  Set work directly and load files.
setwd("/Users/zebinzhang/Desktop/My_Computer/UPPSALA/Capsella_Project")
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

Parents <- c((Accession$populations=="CG" | Accession$populations=="CR"  | Accession$populations=="CO") &
               (Accession$tissue=="flower" | Accession$tissue=="leaf" | Accession$tissue=="root"))

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
### Calculate the library sizes 
# Library size for Cbp = (Co + Cg)*0.5:
libColsums <- colSums(RNAreadsCount + 1, na.rm = T) # add one because later we will have to log values.
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
### Count per million (CPM) normilazation
CPMnorm <- RNAreadsCount
for (i in c(1:c(length(data.lib.size)))){
  CPMnorm[,i] <- ((RNAreadsCount[,i])/data.lib.size[i])*1000000 # why needs to +1 in here?
}
head(CPMnorm)
dim(CPMnorm)

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

row.names(CPMnorm.Pop) <- as.character(CPMnorm.Pop$Gene)
CPMnorm.Pop$Gene <- NULL

head(CPMnorm.Pop)
dim(CPMnorm.Pop)
dim(na.omit(CPMnorm.Pop))
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

###############################
## PCA in Flowers
### Define flower
head(CPMnorm.Pop)
Flower <- c("CR_F","CG_F","CO_F")
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
pca.flower.data$population <- c(c("CR","CG","CO"))
pca.flower.data$population <- factor(pca.flower.data$population, levels = c("CG","CR","CO"))
pca.flower.data$species <- c(c("CR","CG","CO"))


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
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("indianred1","darkviolet","blue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.flower.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.flower.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Flowers")

###############################
## PCA in Leaves
### Define leaf
head(CPMnorm.Pop)
Leaf <- c("CR_L","CG_L","CO_L")
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
pca.leaf.data$population <- c(c("CR","CG","CO"))
pca.leaf.data$population <- factor(pca.leaf.data$population, levels = c("CG","CR","CO"))
pca.leaf.data$species <- c(c("CR","CG","CO"))


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
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("indianred1","darkviolet","blue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.leaf.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.leaf.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Leaves")
###############################
## PCA in Roots
### Define root
head(CPMnorm.Pop)
Root <- c("CR_R","CG_R","CO_R")
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
pca.root.data$population <- c(c("CR","CG","CO"))
pca.root.data$population <- factor(pca.root.data$population, levels = c("CG","CR","CO"))
pca.root.data$species <- c(c("CR","CG","CO"))


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
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("indianred1","darkviolet","blue"))  +
  scale_fill_manual(values = c("indianred1","darkviolet","blue")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Gene expression in Roots")

