setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

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
Accession <- read.csv("InputData/28genomes_annot_R.csv")
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
######################################################
#### Load phenotypic data
Diploids.Pheno <- read.table("InputData/DiploidsPhenotypicFile.txt", header = T, sep = "\t", stringsAsFactors = F)
head(Diploids.Pheno)
Diploids.Pheno$species <- factor(Diploids.Pheno$species, levels = c("CG","CR","CO"))
str(Diploids.Pheno)
## Flower
F.cgco <- Diploids.Pheno[((Diploids.Pheno$tissue == "flower") & 
                            ((Diploids.Pheno$species == "CG") | (Diploids.Pheno$species == "CO"))),]
row.names(F.cgco) <- F.cgco$NewName
str(F.cgco)
F.cgco
F.cgcr <- rbind(Diploids.Pheno[(Diploids.Pheno$tissue == "flower") & (Diploids.Pheno$species == "CG"),],
                Diploids.Pheno[(Diploids.Pheno$tissue == "flower") & (Diploids.Pheno$species == "CR"),])
row.names(F.cgcr) <- F.cgcr$NewName
F.cgcr
str(F.cgcr)
F.crco <- Diploids.Pheno[((Diploids.Pheno$tissue == "flower") & 
                            ((Diploids.Pheno$species == "CR") | (Diploids.Pheno$species == "CO"))),]
row.names(F.crco) <- F.crco$NewName
F.crco
str(F.crco)

## Leaf
L.cgco <- Diploids.Pheno[((Diploids.Pheno$tissue == "leaf") & 
                            ((Diploids.Pheno$species == "CG") | (Diploids.Pheno$species == "CO"))),]
row.names(L.cgco) <- L.cgco$NewName
str(L.cgco)
L.cgcr <- rbind(Diploids.Pheno[(Diploids.Pheno$tissue == "leaf") & (Diploids.Pheno$species == "CG"),],
                Diploids.Pheno[(Diploids.Pheno$tissue == "leaf") & (Diploids.Pheno$species == "CR"),])
row.names(L.cgcr) <- L.cgcr$NewName
str(L.cgcr)
L.crco <- Diploids.Pheno[((Diploids.Pheno$tissue == "leaf") & ((Diploids.Pheno$species == "CR") | (Diploids.Pheno$species == "CO"))),]
row.names(L.crco) <- L.crco$NewName
str(L.crco)

## Root
R.cgco <- Diploids.Pheno[((Diploids.Pheno$tissue == "root") & ((Diploids.Pheno$species == "CG") | (Diploids.Pheno$species == "CO"))),]
row.names(R.cgco) <- R.cgco$NewName
str(R.cgco)
R.cgcr <- rbind(Diploids.Pheno[(Diploids.Pheno$tissue == "root") & (Diploids.Pheno$species == "CG"),],
                Diploids.Pheno[(Diploids.Pheno$tissue == "root") & (Diploids.Pheno$species == "CR"),])
row.names(R.cgcr) <- R.cgcr$NewName
str(R.cgcr)
R.crco <- Diploids.Pheno[((Diploids.Pheno$tissue == "root") & ((Diploids.Pheno$species == "CR") | (Diploids.Pheno$species == "CO"))),]
row.names(R.crco) <- R.crco$NewName
str(R.crco)

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
# Run following commands if analysis include Cbp.                     
#                 (rowSums(raw.RC[,ASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Cg_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Cg_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Cg_F] >= 1, na.rm = T) >= geneMis) &
#             
#                 (rowSums(raw.RC[,ASI_Co_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Co_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Co_F] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Co_F] >= 1, na.rm = T) >= geneMis) &
#             
#                 (rowSums(raw.RC[,ASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Cg_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Cg_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Cg_L] >= 1, na.rm = T) >= geneMis) &
#             
#                 (rowSums(raw.RC[,ASI_Co_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Co_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Co_L] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Co_L] >= 1, na.rm = T) >= geneMis) &
#             
#                 (rowSums(raw.RC[,ASI_Cg_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Cg_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Cg_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Cg_R] >= 1, na.rm = T) >= geneMis) &
#                 
#                 (rowSums(raw.RC[,ASI_Co_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,EUR_Co_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,ME_Co_R] >= 1, na.rm = T) >= geneMis) &
#                 (rowSums(raw.RC[,CASI_Co_R] >= 1, na.rm = T) >= geneMis)
)

raw.RC.Keep <- raw.RC[geneKeep,]
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.Keep)[1]))

### DE with EdgeR
# Create EdgeR dataset
# dge = differentical gene expression
F.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_F | CO_F)],         group = F.cgco$species)
F.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_F], raw.RC.Keep[,CR_F]), group = F.cgcr$species)
F.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_F | CO_F)],         group = F.crco$species)

L.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_L | CO_L)],         group = L.cgco$species)
L.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_L], raw.RC.Keep[,CR_L]), group = L.cgcr$species)
L.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_L | CO_L)],         group = L.crco$species)

R.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_R | CO_R)],         group = R.cgco$species)
R.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_R], raw.RC.Keep[,CR_R]), group = R.cgcr$species)
R.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_R | CO_R)],         group = R.crco$species)
# Normalize by total count 
# TMM normalization
F.cgco.dge <- calcNormFactors(F.cgco.dge, method = "TMM")
F.cgcr.dge <- calcNormFactors(F.cgcr.dge, method = "TMM")
F.crco.dge <- calcNormFactors(F.crco.dge, method = "TMM")

L.cgco.dge <- calcNormFactors(L.cgco.dge, method = "TMM")
L.cgcr.dge <- calcNormFactors(L.cgcr.dge, method = "TMM")
L.crco.dge <- calcNormFactors(L.crco.dge, method = "TMM")

R.cgco.dge <- calcNormFactors(R.cgco.dge, method = "TMM")
R.cgcr.dge <- calcNormFactors(R.cgcr.dge, method = "TMM")
R.crco.dge <- calcNormFactors(R.crco.dge, method = "TMM")
# Generate the barplot of library size
barplot(F.cgco.dge$samples$lib.size, names.arg=F.cgco$NewName,
        las =2, border = F, ylab = "Library size (million reads) in Flower", cex.names=0.5)
# Data Exploration
# Before proceeding with the computations for differential expression, it is possible to produce a plot showing 
# the sample relations based on multidimensional scaling. This is something that we will cover in much more
# detail in a later lecture. The basic premise is that we make a plot so samples which are similar are near to
# each other in the plot while samples that are dissimilar are far from each other. 

plotMDS(F.cgco.dge, method="bcv", col=as.numeric(F.cgco.dge$samples$group))
legend("topright", as.character(unique(F.cgco.dge$samples$group)), col=1:3, pch=20)
plotMDS(F.cgcr.dge, method="bcv", col=as.numeric(F.cgcr.dge$samples$group))
plotMDS(F.crco.dge, method="bcv", col=as.numeric(F.crco.dge$samples$group))

plotMDS(L.cgco.dge, method="bcv", col=as.numeric(L.cgco.dge$samples$group))
legend("topright", as.character(unique(L.cgco.dge$samples$group)), col=1:3, pch=20)
plotMDS(L.cgcr.dge, method="bcv", col=as.numeric(L.cgcr.dge$samples$group))
plotMDS(L.crco.dge, method="bcv", col=as.numeric(L.crco.dge$samples$group))

plotMDS(R.cgco.dge, method="bcv", col=as.numeric(R.cgco.dge$samples$group))
legend("topright", as.character(unique(R.cgco.dge$samples$group)), col=1:3, pch=20)
plotMDS(R.cgcr.dge, method="bcv", col=as.numeric(R.cgcr.dge$samples$group))
plotMDS(R.crco.dge, method="bcv", col=as.numeric(R.crco.dge$samples$group))

# Extract TMM values.
TMM.F.cgco <- estimateCommonDisp(F.cgco.dge)

# Create the contrast matrix
F.cgco.dge.mat <- model.matrix(~ 0 + F.cgco.dge$samples$group)
colnames(F.cgco.dge.mat) <- levels(F.cgco.dge$samples$group)
F.cgcr.dge.mat <- model.matrix(~ 0 + F.cgcr.dge$samples$group)
colnames(F.cgcr.dge.mat) <- levels(F.cgcr.dge$samples$group)
F.crco.dge.mat <- model.matrix(~ 0 + F.crco.dge$samples$group)
colnames(F.crco.dge.mat) <- levels(F.crco.dge$samples$group)

L.cgco.dge.mat <- model.matrix(~ 0 + L.cgco.dge$samples$group)
colnames(L.cgco.dge.mat) <- levels(L.cgco.dge$samples$group)
L.cgcr.dge.mat <- model.matrix(~ 0 + L.cgcr.dge$samples$group)
colnames(L.cgcr.dge.mat) <- levels(L.cgcr.dge$samples$group)
L.crco.dge.mat <- model.matrix(~ 0 + L.crco.dge$samples$group)
colnames(L.crco.dge.mat) <- levels(L.crco.dge$samples$group)

R.cgco.dge.mat <- model.matrix(~ 0 + R.cgco.dge$samples$group)
colnames(R.cgco.dge.mat) <- levels(R.cgco.dge$samples$group)
R.cgcr.dge.mat <- model.matrix(~ 0 + R.cgcr.dge$samples$group)
colnames(R.cgcr.dge.mat) <- levels(R.cgcr.dge$samples$group)
R.crco.dge.mat <- model.matrix(~ 0 + R.crco.dge$samples$group)
colnames(R.crco.dge.mat) <- levels(R.crco.dge$samples$group)

# Estimate dispersion parameter for GLM
# Flower
F.cgco.dge <- estimateGLMCommonDisp(F.cgco.dge, F.cgco.dge.mat)
F.cgco.dge <- estimateGLMTrendedDisp(F.cgco.dge, F.cgco.dge.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
# I chose the method="power" after looking at the methods that are offered: bin.spline (default if number of 
# tags is > 200), power (default otherwise), bin.loess, and spline. We have 17,396 tags, so the default is 
# bin.spline. When I used the bin.spline method, it was no better than estimating a common dispersion so I 
# instead used power.
F.cgco.dge <- estimateGLMTagwiseDisp(F.cgco.dge, F.cgco.dge.mat)

F.cgcr.dge <- estimateGLMCommonDisp(F.cgcr.dge, F.cgcr.dge.mat)
F.cgcr.dge <- estimateGLMTrendedDisp(F.cgcr.dge, F.cgcr.dge.mat, method="power")
F.cgcr.dge <- estimateGLMTagwiseDisp(F.cgcr.dge, F.cgcr.dge.mat)

F.crco.dge <- estimateGLMCommonDisp(F.crco.dge, F.crco.dge.mat)
F.crco.dge <- estimateGLMTrendedDisp(F.crco.dge, F.crco.dge.mat, method="power")
F.crco.dge <- estimateGLMTagwiseDisp(F.crco.dge, F.crco.dge.mat)

# Leaf
L.cgco.dge <- estimateGLMCommonDisp(L.cgco.dge, L.cgco.dge.mat)
L.cgco.dge <- estimateGLMTrendedDisp(L.cgco.dge, L.cgco.dge.mat, method="power")
L.cgco.dge <- estimateGLMTagwiseDisp(L.cgco.dge, L.cgco.dge.mat)

L.cgcr.dge <- estimateGLMCommonDisp(L.cgcr.dge, L.cgcr.dge.mat)
L.cgcr.dge <- estimateGLMTrendedDisp(L.cgcr.dge, L.cgcr.dge.mat, method="power")
L.cgcr.dge <- estimateGLMTagwiseDisp(L.cgcr.dge, L.cgcr.dge.mat)

L.crco.dge <- estimateGLMCommonDisp(L.crco.dge, L.crco.dge.mat)
L.crco.dge <- estimateGLMTrendedDisp(L.crco.dge, L.crco.dge.mat, method="power")
L.crco.dge <- estimateGLMTagwiseDisp(L.crco.dge, L.crco.dge.mat)

# Root
R.cgco.dge <- estimateGLMCommonDisp(R.cgco.dge, R.cgco.dge.mat)
R.cgco.dge <- estimateGLMTrendedDisp(R.cgco.dge, R.cgco.dge.mat, method="power")
R.cgco.dge <- estimateGLMTagwiseDisp(R.cgco.dge, R.cgco.dge.mat)

R.cgcr.dge <- estimateGLMCommonDisp(R.cgcr.dge, R.cgcr.dge.mat)
R.cgcr.dge <- estimateGLMTrendedDisp(R.cgcr.dge, R.cgcr.dge.mat, method="power")
R.cgcr.dge <- estimateGLMTagwiseDisp(R.cgcr.dge, R.cgcr.dge.mat)

R.crco.dge <- estimateGLMCommonDisp(R.crco.dge, R.crco.dge.mat)
R.crco.dge <- estimateGLMTrendedDisp(R.crco.dge, R.crco.dge.mat, method="power")
R.crco.dge <- estimateGLMTagwiseDisp(R.crco.dge, R.crco.dge.mat)


# Plot mean-variance
plotBCV(F.cgco.dge)
plotBCV(F.cgcr.dge)
plotBCV(F.crco.dge)
plotBCV(L.cgco.dge)
plotBCV(L.cgcr.dge)
plotBCV(L.crco.dge)
plotBCV(R.cgco.dge)
plotBCV(R.cgcr.dge)
plotBCV(R.crco.dge)
##########
## DEG analysis 
# Set threshold
p.threshold <- 0.01
## edgeR ##
# Model fitting
F.cgco.fit <- glmFit(F.cgco.dge, F.cgco.dge.mat)
F.cgcr.fit <- glmFit(F.cgcr.dge, F.cgcr.dge.mat)
F.crco.fit <- glmFit(F.crco.dge, F.crco.dge.mat)

L.cgco.fit <- glmFit(L.cgco.dge, L.cgco.dge.mat)
L.cgcr.fit <- glmFit(L.cgcr.dge, L.cgcr.dge.mat)
L.crco.fit <- glmFit(L.crco.dge, L.crco.dge.mat)

R.cgco.fit <- glmFit(R.cgco.dge, R.cgco.dge.mat)
R.cgcr.fit <- glmFit(R.cgcr.dge, R.cgcr.dge.mat)
R.crco.fit <- glmFit(R.crco.dge, R.crco.dge.mat)

# Differential expression
# Flower
F.cgco.contrast <- makeContrasts(CG - CO, levels = F.cgco.dge.mat)
F.cgco.lrt <- glmLRT(F.cgco.fit, contrast = F.cgco.contrast)

F.cgcr.contrast <- makeContrasts(CG - CR, levels = F.cgcr.dge.mat)
F.cgcr.lrt <- glmLRT(F.cgcr.fit, contrast = F.cgcr.contrast)

F.crco.contrast <- makeContrasts(CR - CO, levels = F.crco.dge.mat)
F.crco.lrt <- glmLRT(F.crco.fit, contrast = F.crco.contrast)
# Leaf
L.cgco.contrast <- makeContrasts(CG - CO, levels = L.cgco.dge.mat)
L.cgco.lrt <- glmLRT(L.cgco.fit, contrast = L.cgco.contrast)

L.cgcr.contrast <- makeContrasts(CG - CR, levels = L.cgcr.dge.mat)
L.cgcr.lrt <- glmLRT(L.cgcr.fit, contrast = L.cgcr.contrast)

L.crco.contrast <- makeContrasts(CR - CO, levels = L.crco.dge.mat)
L.crco.lrt <- glmLRT(L.crco.fit, contrast = L.crco.contrast)
# Root
R.cgco.contrast <- makeContrasts(CG - CO, levels = R.cgco.dge.mat)
R.cgco.lrt <- glmLRT(R.cgco.fit, contrast = R.cgco.contrast)

R.cgcr.contrast <- makeContrasts(CG - CR, levels = R.cgcr.dge.mat)
R.cgcr.lrt <- glmLRT(R.cgcr.fit, contrast = R.cgcr.contrast)

R.crco.contrast <- makeContrasts(CR - CO, levels = R.crco.dge.mat)
R.crco.lrt <- glmLRT(R.crco.fit, contrast = R.crco.contrast)

# Access results tables
F.cgco.table <- F.cgco.lrt$table
F.cgcr.table <- F.cgcr.lrt$table
F.crco.table <- F.crco.lrt$table

L.cgco.table <- L.cgco.lrt$table
L.cgcr.table <- L.cgcr.lrt$table
L.crco.table <- L.crco.lrt$table

R.cgco.table <- R.cgco.lrt$table
R.cgcr.table <- R.cgcr.lrt$table
R.crco.table <- R.crco.lrt$table

# Toptags Result with FDR
F.cgco.topTags <- topTags(F.cgco.lrt, n = nrow(F.cgco.lrt$table))$table
F.cgcr.topTags <- topTags(F.cgcr.lrt, n = nrow(F.cgcr.lrt$table))$table
F.crco.topTags <- topTags(F.crco.lrt, n = nrow(F.crco.lrt$table))$table

L.cgco.topTags <- topTags(L.cgco.lrt, n = nrow(L.cgco.lrt$table))$table
L.cgcr.topTags <- topTags(L.cgcr.lrt, n = nrow(L.cgcr.lrt$table))$table
L.crco.topTags <- topTags(L.crco.lrt, n = nrow(L.crco.lrt$table))$table

R.cgco.topTags <- topTags(R.cgco.lrt, n = nrow(R.cgco.lrt$table))$table
R.cgcr.topTags <- topTags(R.cgcr.lrt, n = nrow(R.cgcr.lrt$table))$table
R.crco.topTags <- topTags(R.crco.lrt, n = nrow(R.crco.lrt$table))$table

# Extract the DE genes with FDR < 0.01
# FC > 1 & FDR ≤ 0.01
# DE.F.cgco <- F.cgco.topTags[(abs(F.cgco.topTags$logFC) > 1) & (F.cgco.topTags$FDR <= 0.01),]
# DE.F.cgcr <- F.cgcr.topTags[(abs(F.cgcr.topTags$logFC) > 1) & (F.cgcr.topTags$FDR <= 0.01),]
# DE.F.crco <- F.crco.topTags[(abs(F.crco.topTags$logFC) > 1) & (F.crco.topTags$FDR <= 0.01),]
# 
# DE.L.cgco <- L.cgco.topTags[(abs(L.cgco.topTags$logFC) > 1) & (L.cgco.topTags$FDR <= 0.01),]
# DE.L.cgcr <- L.cgcr.topTags[(abs(L.cgcr.topTags$logFC) > 1) & (L.cgcr.topTags$FDR <= 0.01),]
# DE.L.crco <- L.crco.topTags[(abs(L.crco.topTags$logFC) > 1) & (L.crco.topTags$FDR <= 0.01),]
# 
# DE.R.cgco <- R.cgco.topTags[(abs(R.cgco.topTags$logFC) > 1) & (R.cgco.topTags$FDR <= 0.01),]
# DE.R.cgcr <- R.cgcr.topTags[(abs(R.cgcr.topTags$logFC) > 1) & (R.cgcr.topTags$FDR <= 0.01),]
# DE.R.crco <- R.crco.topTags[(abs(R.crco.topTags$logFC) > 1) & (R.crco.topTags$FDR <= 0.01),]

# FDR ≤ 0.01
DE.F.cgco <- F.cgco.topTags[(F.cgco.topTags$FDR <= 0.01),]
DE.F.cgcr <- F.cgcr.topTags[(F.cgcr.topTags$FDR <= 0.01),]
DE.F.crco <- F.crco.topTags[(F.crco.topTags$FDR <= 0.01),]

DE.L.cgco <- L.cgco.topTags[(L.cgco.topTags$FDR <= 0.01),]
DE.L.cgcr <- L.cgcr.topTags[(L.cgcr.topTags$FDR <= 0.01),]
DE.L.crco <- L.crco.topTags[(L.crco.topTags$FDR <= 0.01),]

DE.R.cgco <- R.cgco.topTags[(R.cgco.topTags$FDR <= 0.01),]
DE.R.cgcr <- R.cgcr.topTags[(R.cgcr.topTags$FDR <= 0.01),]
DE.R.crco <- R.crco.topTags[(R.crco.topTags$FDR <= 0.01),]

message("Number of DE genes between Cg and Co in Flower ", dim(DE.F.cgco)[1])
message("Number of DE genes between Cg and Cr in Flower ", dim(DE.F.cgcr)[1])
message("Number of DE genes between Cr and Co in Flower ", dim(DE.F.crco)[1])

message("Number of DE genes between Cg and Co in Leaf ", dim(DE.L.cgco)[1])
message("Number of DE genes between Cg and Cr in Leaf ", dim(DE.L.cgcr)[1])
message("Number of DE genes between Cr and Co in Leaf ", dim(DE.L.crco)[1])

message("Number of DE genes between Cg and Co in Root ", dim(DE.R.cgco)[1])
message("Number of DE genes between Cg and Cr in Root ", dim(DE.R.cgcr)[1])
message("Number of DE genes between Cr and Co in Root ", dim(DE.R.crco)[1])

# write.table(DE.F.cgco, file = "OutputData/DE_CgCo_F genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.F.cgcr, file = "OutputData/DE_CgCr_F genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.F.crco, file = "OutputData/DE_CrCo_F genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# 
# write.table(DE.L.cgco, file = "OutputData/DE_CgCo_L genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.L.cgcr, file = "OutputData/DE_CgCr_L genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.L.crco, file = "OutputData/DE_CrCo_L genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# 
# write.table(DE.R.cgco, file = "OutputData/DE_CgCo_R genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.R.cgcr, file = "OutputData/DE_CgCr_R genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)
# write.table(DE.R.crco, file = "OutputData/DE_CrCo_R genes FDR 0.01.txt", append = F, quote = F, sep = "\t", 
#             row.names = T)


############
# # BH "Benjamini-Hochberg"
# Flower.sig.edgeR <- decideTestsDGE(Flower.lrt.edgeR, adjust.method="BH", p.value = p.threshold)
# # add BH filter result to Access results tables.
# Flower.topTags.edgeR$BH <- 0
# Flower.topTags.edgeR[which(Flower.sig.edgeR != 0),]$BH <- 1
# # Calculate how many genes passed the BH filter, result is 2,134
# dim(Flower.topTags.edgeR[Flower.topTags.edgeR$BH ==1,])
# # Sort the logFC, and pick the first top 5 percentage value, is 5.2
# sort(abs(Flower.topTags.edgeR$logFC), decreasing = TRUE)
# Flower.topTags.edgeR$Gene <- row.names(Flower.topTags.edgeR)
# head(Flower.topTags.edgeR)
# dim(Flower.topTags.edgeR)
# dim(Flower.topTags.edgeR[((Flower.topTags.edgeR$FDR < 0.05) & (abs(Flower.topTags.edgeR$logFC) > 5.2)),])
# FlowerKeep.EdgeR <- Flower.topTags.edgeR[((Flower.topTags.edgeR$FDR < 0.05) & (abs(Flower.topTags.edgeR$logFC) > 5.2)),]$Gene
############

##### Make a basic volcano plot #####
par()              # view current settings
opar <- par()
par(mfrow=c(3,3))
# F.cgco
with(F.cgco.topTags, plot(logFC, 
               -log10(FDR), 
               pch=20, 
               main="Gene expression of Cg vs. Co in Flower"))

with(subset(F.cgco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="firebrick1"))
# F.cgcr
with(F.cgcr.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cg vs. Cr in Flower"))

with(subset(F.cgcr.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="forestgreen"))
# F.crco
with(F.crco.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cr vs. Co in Flower"))

with(subset(F.crco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="purple4"))
# L.cgco
with(L.cgco.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cg vs. Co in Leaf"))

with(subset(L.cgco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="firebrick1"))
# L.cgcr
with(L.cgcr.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cg vs. Cr in Leaf"))

with(subset(L.cgcr.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="forestgreen"))
# L.crco
with(L.crco.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cr vs. Co in Leaf"))

with(subset(L.crco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="purple4"))
# R.cgco
with(R.cgco.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cg vs. Co in Root"))

with(subset(R.cgco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="firebrick1"))
# R.cgcr
with(R.cgcr.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cg vs. Cr in Root"))

with(subset(R.cgcr.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="forestgreen"))
# R.crco
with(R.crco.topTags, plot(logFC, 
                          -log10(FDR), 
                          pch=20, 
                          main="Gene expression of Cr vs. Co in Root"))

with(subset(R.crco.topTags,  abs(logFC) > 1 & FDR <= 0.01), 
     points(logFC, 
            -log10(FDR), 
            pch=20, 
            col="purple4"))
par(opar)
dev.off()

#### Histogram ####
par()              # view current settings
opar <- par()
par(mfrow=c(3,3))

hist(F.cgco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Co in flower", xlab="logFC", ylab = "Counts")
hist(F.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Cr in flower", xlab="logFC", ylab = "Counts")
hist(F.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cr vs. Co in flower", xlab="logFC", ylab = "Counts")

hist(L.cgco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Co in leaf", xlab="logFC", ylab = "Counts")
hist(L.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Cr in leaf", xlab="logFC", ylab = "Counts")
hist(L.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cr vs. Co in leaf", xlab="logFC", ylab = "Counts")

hist(R.cgco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Co in root", xlab="logFC", ylab = "Counts")
hist(R.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cg vs. Cr in root", xlab="logFC", ylab = "Counts")
hist(R.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes of Cr vs. Co in root", xlab="logFC", ylab = "Counts")

par(opar)
dev.off()

#### ggscatter ####
library(ggpubr)
#install.packages("ggExtra")
library(ggExtra)

F.cgco.topTags$logFDR <- -log10(F.cgco.topTags$FDR)
F.cgco.topTags$Tissue <- "Flower"
F.cgco.topTags$DE <- "Not DE"
F.cgco.topTags[(F.cgco.topTags$logFC > 1) & (F.cgco.topTags$FDR <= 0.01),]$DE <- "Up"
F.cgco.topTags[(F.cgco.topTags$logFC < -1) & (F.cgco.topTags$FDR <= 0.01),]$DE <- "Down"
str(F.cgco.topTags)
F.cgco.topTags$DE <- factor(F.cgco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(F.cgco.topTags)


F.cgcr.topTags$logFDR <- -log10(F.cgcr.topTags$FDR)
F.cgcr.topTags$Tissue <- "Flower"
F.cgcr.topTags$DE <- "Not DE"
F.cgcr.topTags[(F.cgcr.topTags$logFC > 1) & (F.cgcr.topTags$FDR <= 0.01),]$DE <- "Up"
F.cgcr.topTags[(F.cgcr.topTags$logFC < -1) & (F.cgcr.topTags$FDR <= 0.01),]$DE <- "Down"
str(F.cgcr.topTags)
F.cgcr.topTags$DE <- factor(F.cgcr.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(F.cgcr.topTags)

F.crco.topTags$logFDR <- -log10(F.crco.topTags$FDR)
F.crco.topTags$Tissue <- "Flower"
F.crco.topTags$DE <- "Not DE"
F.crco.topTags[(F.crco.topTags$logFC > 1) & (F.crco.topTags$FDR <= 0.01),]$DE <- "Up"
F.crco.topTags[(F.crco.topTags$logFC < -1) & (F.crco.topTags$FDR <= 0.01),]$DE <- "Down"
str(F.crco.topTags)
F.crco.topTags$DE <- factor(F.crco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(F.crco.topTags)

L.cgco.topTags$logFDR <- -log10(L.cgco.topTags$FDR)
L.cgco.topTags$Tissue <- "Leaf"
L.cgco.topTags$DE <- "Not DE"
L.cgco.topTags[(L.cgco.topTags$logFC > 1) & (L.cgco.topTags$FDR <= 0.01),]$DE <- "Up"
L.cgco.topTags[(L.cgco.topTags$logFC < -1) & (L.cgco.topTags$FDR <= 0.01),]$DE <- "Down"
str(L.cgco.topTags)
L.cgco.topTags$DE <- factor(L.cgco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(L.cgco.topTags)

L.cgcr.topTags$logFDR <- -log10(L.cgcr.topTags$FDR)
L.cgcr.topTags$Tissue <- "Leaf"
L.cgcr.topTags$DE <- "Not DE"
L.cgcr.topTags[(L.cgcr.topTags$logFC > 1) & (L.cgcr.topTags$FDR <= 0.01),]$DE <- "Up"
L.cgcr.topTags[(L.cgcr.topTags$logFC < -1) & (L.cgcr.topTags$FDR <= 0.01),]$DE <- "Down"
str(L.cgcr.topTags)
L.cgcr.topTags$DE <- factor(L.cgcr.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(L.cgcr.topTags)

L.crco.topTags$logFDR <- -log10(L.crco.topTags$FDR)
L.crco.topTags$Tissue <- "Leaf"
L.crco.topTags$DE <- "Not DE"
L.crco.topTags[(L.crco.topTags$logFC > 1) & (L.crco.topTags$FDR <= 0.01),]$DE <- "Up"
L.crco.topTags[(L.crco.topTags$logFC < -1) & (L.crco.topTags$FDR <= 0.01),]$DE <- "Down"
str(L.crco.topTags)
L.crco.topTags$DE <- factor(L.crco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(L.crco.topTags)

R.cgco.topTags$logFDR <- -log10(R.cgco.topTags$FDR)
R.cgco.topTags$Tissue <- "Root"
R.cgco.topTags$DE <- "Not DE"
R.cgco.topTags[(R.cgco.topTags$logFC > 1) & (R.cgco.topTags$FDR <= 0.01),]$DE <- "Up"
R.cgco.topTags[(R.cgco.topTags$logFC < -1) & (R.cgco.topTags$FDR <= 0.01),]$DE <- "Down"
str(R.cgco.topTags)
R.cgco.topTags$DE <- factor(R.cgco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(R.cgco.topTags)

R.cgcr.topTags$logFDR <- -log10(R.cgcr.topTags$FDR)
R.cgcr.topTags$Tissue <- "Root"
R.cgcr.topTags$DE <- "Not DE"
R.cgcr.topTags[(R.cgcr.topTags$logFC > 1) & (R.cgcr.topTags$FDR <= 0.01),]$DE <- "Up"
R.cgcr.topTags[(R.cgcr.topTags$logFC < -1) & (R.cgcr.topTags$FDR <= 0.01),]$DE <- "Down"
str(R.cgcr.topTags)
R.cgcr.topTags$DE <- factor(R.cgcr.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(R.cgcr.topTags)

R.crco.topTags$logFDR <- -log10(R.crco.topTags$FDR)
R.crco.topTags$Tissue <- "Root"
R.crco.topTags$DE <- "Not DE"
R.crco.topTags[(R.crco.topTags$logFC > 1) & (R.crco.topTags$FDR <= 0.01),]$DE <- "Up"
R.crco.topTags[(R.crco.topTags$logFC < -1) & (R.crco.topTags$FDR <= 0.01),]$DE <- "Down"
str(R.crco.topTags)
R.crco.topTags$DE <- factor(R.crco.topTags$DE, levels = c("Up", "Down", "Not DE"))
head(R.crco.topTags)

### Message out 
median(F.cgco.topTags[F.cgco.topTags$DE == "Up",]$logFC)
median(F.cgco.topTags[F.cgco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Co in Flower is: ", dim(F.cgco.topTags[F.cgco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Co in Flower is: ", dim(F.cgco.topTags[F.cgco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Co in Flower is: ", (dim(F.cgco.topTags[F.cgco.topTags$DE == "Up",])[1] / dim(F.cgco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Co in Flower is: ", (dim(F.cgco.topTags[F.cgco.topTags$DE == "Down",])[1] / dim(F.cgco.topTags)[1])*100)

median(F.cgcr.topTags[F.cgcr.topTags$DE == "Up",]$logFC)
median(F.cgcr.topTags[F.cgcr.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Cr in Flower is: ", dim(F.cgcr.topTags[F.cgcr.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Cr in Flower is: ", dim(F.cgcr.topTags[F.cgcr.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Cr in Flower is: ", (dim(F.cgcr.topTags[F.cgcr.topTags$DE == "Up",])[1] / dim(F.cgcr.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Cr in Flower is: ", (dim(F.cgcr.topTags[F.cgcr.topTags$DE == "Down",])[1] / dim(F.cgcr.topTags)[1])*100)

median(F.crco.topTags[F.crco.topTags$DE == "Up",]$logFC)
median(F.crco.topTags[F.crco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cr vs. Co in Flower is: ", dim(F.crco.topTags[F.crco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cr vs. Co in Flower is: ", dim(F.crco.topTags[F.crco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cr vs. Co in Flower is: ", (dim(F.crco.topTags[F.crco.topTags$DE == "Up",])[1] / dim(F.crco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cr vs. Co in Flower is: ", (dim(F.crco.topTags[F.crco.topTags$DE == "Down",])[1] / dim(F.crco.topTags)[1])*100)

median(L.cgco.topTags[L.cgco.topTags$DE == "Up",]$logFC)
median(L.cgco.topTags[L.cgco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Co in Leaf is: ", dim(L.cgco.topTags[L.cgco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Co in Leaf is: ", dim(L.cgco.topTags[L.cgco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Co in Leaf is: ", (dim(L.cgco.topTags[L.cgco.topTags$DE == "Up",])[1] / dim(L.cgco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Co in Leaf is: ", (dim(L.cgco.topTags[L.cgco.topTags$DE == "Down",])[1] / dim(L.cgco.topTags)[1])*100)

median(L.cgcr.topTags[L.cgcr.topTags$DE == "Up",]$logFC)
median(L.cgcr.topTags[L.cgcr.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Cr in Leaf is: ", dim(L.cgcr.topTags[L.cgcr.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Cr in Leaf is: ", dim(L.cgcr.topTags[L.cgcr.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Cr in Leaf is: ", (dim(L.cgcr.topTags[L.cgcr.topTags$DE == "Up",])[1] / dim(L.cgcr.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Cr in Leaf is: ", (dim(L.cgcr.topTags[L.cgcr.topTags$DE == "Down",])[1] / dim(L.cgcr.topTags)[1])*100)

median(L.crco.topTags[L.crco.topTags$DE == "Up",]$logFC)
median(L.crco.topTags[L.crco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cr vs. Co in Leaf is: ", dim(L.crco.topTags[L.crco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cr vs. Co in Leaf is: ", dim(L.crco.topTags[L.crco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cr vs. Co in Leaf is: ", (dim(L.crco.topTags[L.crco.topTags$DE == "Up",])[1] / dim(L.crco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cr vs. Co in Leaf is: ", (dim(L.crco.topTags[L.crco.topTags$DE == "Down",])[1] / dim(L.crco.topTags)[1])*100)

median(R.cgco.topTags[R.cgco.topTags$DE == "Up",]$logFC)
median(R.cgco.topTags[R.cgco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Co in Root is: ", dim(R.cgco.topTags[R.cgco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Co in Root is: ", dim(R.cgco.topTags[R.cgco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Co in Root is: ", (dim(R.cgco.topTags[R.cgco.topTags$DE == "Up",])[1] / dim(R.cgco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Co in Root is: ", (dim(R.cgco.topTags[R.cgco.topTags$DE == "Down",])[1] / dim(R.cgco.topTags)[1])*100)

median(R.cgcr.topTags[R.cgcr.topTags$DE == "Up",]$logFC)
median(R.cgcr.topTags[R.cgcr.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cg vs. Cr in Root is: ", dim(R.cgcr.topTags[R.cgcr.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cg vs. Cr in Root is: ", dim(R.cgcr.topTags[R.cgcr.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cg vs. Cr in Root is: ", (dim(R.cgcr.topTags[R.cgcr.topTags$DE == "Up",])[1] / dim(R.cgcr.topTags)[1])*100)
message("Percentage of Down regulated genes of Cg vs. Cr in Root is: ", (dim(R.cgcr.topTags[R.cgcr.topTags$DE == "Down",])[1] / dim(R.cgcr.topTags)[1])*100)

median(R.crco.topTags[R.crco.topTags$DE == "Up",]$logFC)
median(R.crco.topTags[R.crco.topTags$DE == "Down",]$logFC)
message("Numbers of Up regulated genes of Cr vs. Co in Root is: ", dim(R.crco.topTags[R.crco.topTags$DE == "Up",])[1])
message("Numbers of Down regulated genes of Cr vs. Co in Root is: ", dim(R.crco.topTags[R.crco.topTags$DE == "Down",])[1])
message("Percentage of Up regulated genes of Cr vs. Co in Root is: ", (dim(R.crco.topTags[R.crco.topTags$DE == "Up",])[1] / dim(R.crco.topTags)[1])*100)
message("Percentage of Down regulated genes of Cr vs. Co in Root is: ", (dim(R.crco.topTags[R.crco.topTags$DE == "Down",])[1] / dim(R.crco.topTags)[1])*100)


############## 
## Scatter plot
## p <- ggscatter(F.cgco.topTags, x = "logFC", y = "logPValue",
##           color = "DE", palette = "jco",         # Color by groups "DE"
##           size = 1, alpha = 0.6      # change the plot size and transparency
## 
## ) 
## # Marginal plot  
## xplot <- ggdensity(F.cgco.topTags, "logFC", fill = "DE",
##                    palette = "jco") + clean_theme() + rremove("legend")
## # Cleaning the plots
## p <- p + rremove("legend")
## xplot <- xplot + clean_theme() + rremove("legend")
## # Arranging the plot using cowplot
## library(cowplot)
## plot_grid(xplot,  p, 
##           nrow =  2, # Number of rows in the plot grid.
##           align = "hv", # "hv": align the grid in both directions, "h" horizontally and "v" vertically.
##           rel_widths = c(1, 1), rel_heights = c(1, 2))
## 
## # The problem with the above plots, is the presence of extra spaces between the main plot and the marginal density plots. 
########
### Second pattern
library(ggplot2)
library(cowplot) 
# 1. F.cgco
# Main plot
F_cgco_main <- ggscatter(F.cgco.topTags, x = "logFC", y = "logFDR",
               color = "DE",       # Color by groups "DE"
               palette = c("firebrick1", "sienna1", "gray"), # reset the colour        
               size = 0.2, alpha = 0.5,      # change the plot size and transparency
               title = "Cg vs. Co in Flower", # Add title
               ggtheme = theme_bw() # change background to black pattern
               ) + rremove("legend") + # remove figure legend
               scale_x_continuous(limits = c(-15, 15),  # set the x axis limit
                                  breaks = get_breaks(by = 5, from = -15)) # set the breaks in x axis
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
F_cgco_xdens <- axis_canvas(F_cgco_main, axis = "x") + # marginal will add to main plot along x axis.
                geom_density(data = F.cgco.topTags, aes(x = logFC, fill = DE),
                alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("firebrick1", "sienna1", "gray")) + rremove("legend")

F_cgco_xboxs <- ggboxplot(F.cgco.topTags, x = "DE", y = "logFC",
                color = "DE", fill = "DE", palette = c("firebrick1", "sienna1", "gray"),
                ylim = c(-15,15),
                alpha = 0.5,
                rotate= TRUE,  
                ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") + # add mean value
                rremove("legend")

F_cgco_p1 <- insert_xaxis_grob(F_cgco_main, F_cgco_xdens, grid::unit(.5, "null"), position = "top")
F_cgco_p2 <- insert_xaxis_grob(F_cgco_p1, F_cgco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(F_cgco_p2)

# 2. F.cgcr
# Main plot
F_cgcr_main <- ggscatter(F.cgcr.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("forestgreen", "lightseagreen", "gray"),         
                         size = 0.2, alpha = 0.5,      # change the plot size and transparency
                         title = "Cg vs. Cr in Flower",
                         ggtheme = theme_bw()
                        ) + rremove("legend") +
                        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
F_cgcr_xdens <- axis_canvas(F_cgcr_main, axis = "x")+
                geom_density(data = F.cgcr.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("forestgreen", "lightseagreen", "gray")) + rremove("legend")

F_cgcr_xboxs <- ggboxplot(F.cgcr.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("forestgreen", "lightseagreen", "gray"),
                          alpha = 0.5, 
                          rotate= TRUE,  
                          ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
  stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
  rremove("legend")

F_cgcr_p1 <- insert_xaxis_grob(F_cgcr_main, F_cgcr_xdens, grid::unit(.5, "null"), position = "top")
F_cgcr_p2 <- insert_xaxis_grob(F_cgcr_p1, F_cgcr_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(F_cgcr_p2)

# 3. F.crco
# Main plot
F_crco_main <- ggscatter(F.crco.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("purple4", "orchid1", "gray"),         
                         size = 0.2, alpha = 0.5,      # change the plot size and transparency
                         title = "Cr vs. Co in Flower",
                         ggtheme = theme_bw()
                        ) + rremove("legend") +
                        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
F_crco_xdens <- axis_canvas(F_crco_main, axis = "x")+
                geom_density(data = F.crco.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("purple4", "orchid1", "gray")) + rremove("legend")

F_crco_xboxs <- ggboxplot(F.crco.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("purple4", "orchid1", "gray"),
                          alpha = 0.5, 
                          rotate= TRUE,  
                          ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rremove("legend")

F_crco_p1 <- insert_xaxis_grob(F_crco_main, F_crco_xdens, grid::unit(.5, "null"), position = "top")
F_crco_p2 <- insert_xaxis_grob(F_crco_p1, F_crco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(F_crco_p2)

# 4. L.cgco
# Main plot
L_cgco_main <- ggscatter(L.cgco.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("firebrick1", "sienna1", "gray"),         
                         size = 0.2, alpha = 0.5,      # change the plot size and transparency
                         title = "Cg vs. Co in Leaf",
                         ggtheme = theme_bw()
                        ) + rremove("legend") +
                        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
L_cgco_xdens <- axis_canvas(L_cgco_main, axis = "x")+
                geom_density(data = L.cgco.topTags, aes(x = logFC, fill = DE),
                alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("firebrick1", "sienna1", "gray")) + rremove("legend") 

L_cgco_xboxs <- ggboxplot(L.cgco.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("firebrick1", "sienna1", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

L_cgco_p1 <- insert_xaxis_grob(L_cgco_main, L_cgco_xdens, grid::unit(.5, "null"), position = "top")
L_cgco_p2 <- insert_xaxis_grob(L_cgco_p1, L_cgco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(L_cgco_p2)

# 5. L.cgcr
# Main plot
L_cgcr_main <- ggscatter(L.cgcr.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("forestgreen", "lightseagreen", "gray"),         
                         size = 0.2, alpha = 0.6,      # change the plot size and transparency
                         title = "Cg vs. Cr in Leaf",
                         ggtheme = theme_bw()
                        ) + rremove("legend") +
                        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
L_cgcr_xdens <- axis_canvas(L_cgcr_main, axis = "x")+
                geom_density(data = L.cgcr.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("forestgreen", "lightseagreen", "gray")) + rremove("legend")

L_cgcr_xboxs <- ggboxplot(L.cgcr.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("forestgreen", "lightseagreen", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

L_cgcr_p1 <- insert_xaxis_grob(L_cgcr_main, L_cgcr_xdens, grid::unit(.5, "null"), position = "top")
L_cgcr_p2 <- insert_xaxis_grob(L_cgcr_p1, L_cgcr_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(L_cgcr_p2)

# 6. L.crco
# Main plot
L_crco_main <- ggscatter(L.crco.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("purple4", "orchid1", "gray"),         
                         size = 0.2, alpha = 0.6,      # change the plot size and transparency
                         title = "Cr vs. Co in Leaf",
                         ggtheme = theme_bw()
                        ) + rremove("legend") +
                        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
L_crco_xdens <- axis_canvas(L_crco_main, axis = "x")+
                geom_density(data = L.crco.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("purple4", "orchid1", "gray")) + rremove("legend")

L_crco_xboxs <- ggboxplot(L.crco.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("purple4", "orchid1", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

L_crco_p1 <- insert_xaxis_grob(L_crco_main, L_crco_xdens, grid::unit(.5, "null"), position = "top")
L_crco_p2 <- insert_xaxis_grob(L_crco_p1, L_crco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(L_crco_p2)

# 7. R.cgco
# Main plot
R_cgco_main <- ggscatter(R.cgco.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("firebrick1", "sienna1", "gray"),         
                         size = 0.2, alpha = 0.6,      # change the plot size and transparency
                         title = "Cg vs. Co in Root",
                         ggtheme = theme_bw()
                        ) + 
               rremove("legend") +
               scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
R_cgco_xdens <- axis_canvas(R_cgco_main, axis = "x")+
                geom_density(data = R.cgco.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("firebrick1", "sienna1", "gray")) + rremove("legend")

R_cgco_xboxs <- ggboxplot(R.cgco.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("firebrick1", "sienna1", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

R_cgco_p1 <- insert_xaxis_grob(R_cgco_main, R_cgco_xdens, grid::unit(.5, "null"), position = "top")
R_cgco_p2 <- insert_xaxis_grob(R_cgco_p1, R_cgco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(R_cgco_p2)

# 8. R.cgcr
# Main plot
R_cgcr_main <- ggscatter(R.cgcr.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("forestgreen", "lightseagreen", "gray"),         
                         size = 0.2, alpha = 0.6,      # change the plot size and transparency
                         title = "Cg vs. Cr in Root",
                         ggtheme = theme_bw()
               ) + rremove("legend") +
               scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
R_cgcr_xdens <- axis_canvas(R_cgcr_main, axis = "x")+
                geom_density(data = R.cgcr.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2) +
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("forestgreen", "lightseagreen", "gray")) + rremove("legend")

R_cgcr_xboxs <- ggboxplot(R.cgcr.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("forestgreen", "lightseagreen", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

R_cgcr_p1 <- insert_xaxis_grob(R_cgcr_main, R_cgcr_xdens, grid::unit(.5, "null"), position = "top")
R_cgcr_p2 <- insert_xaxis_grob(R_cgcr_p1, R_cgcr_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(R_cgcr_p2)

# 9. R.crco
# Main plot
R_crco_main <- ggscatter(R.crco.topTags, x = "logFC", y = "logFDR",
                         color = "DE",       # Color by groups "DE"
                         palette = c("purple4", "orchid1", "gray"),         
                         size = 0.2, alpha = 0.6,      # change the plot size and transparency
                         title = "Cr vs. Co in Root",
                         ggtheme = theme_bw()
               ) + rremove("legend") +
               scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15))
# Marginal densities along x axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
R_crco_xdens <- axis_canvas(R_crco_main, axis = "x")+
                geom_density(data = R.crco.topTags, aes(x = logFC, fill = DE),
                             alpha = 0.7, size = 0.2)+
                scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                ggpubr::fill_palette(c("purple4", "orchid1", "gray")) + rremove("legend")

R_crco_xboxs <- ggboxplot(R.crco.topTags, x = "DE", y = "logFC",
                          color = "DE", fill = "DE", palette = c("purple4", "orchid1", "gray"),
                          alpha = 0.5, ggtheme = theme_bw()) +
                scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
                stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
                rotate()  + rremove("legend")

R_crco_p1 <- insert_xaxis_grob(R_crco_main, R_crco_xdens, grid::unit(.5, "null"), position = "top")
R_crco_p2 <- insert_xaxis_grob(R_crco_p1, R_crco_xboxs, grid::unit(.5, "null"), position = "top")
ggdraw(R_crco_p2)

# combine all plot together
ggarrange(F_cgco_p2, F_cgcr_p2, F_crco_p2,
          L_cgco_p2, L_cgcr_p2, L_crco_p2,
          R_cgco_p2, R_cgcr_p2, R_crco_p2,
          labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3, ncol = 3)


# Extract the DE genes between outcrosser and selfer
# Add colume of gene name into dataframe
DE.F.cgco$gene <- rownames(DE.F.cgco)
DE.F.cgcr$gene <- rownames(DE.F.cgcr)
DE.F.crco$gene <- rownames(DE.F.crco)

DE.L.cgco$gene <- rownames(DE.L.cgco)
DE.L.cgcr$gene <- rownames(DE.L.cgcr)
DE.L.crco$gene <- rownames(DE.L.crco)

DE.R.cgco$gene <- rownames(DE.R.cgco)
DE.R.cgcr$gene <- rownames(DE.R.cgcr)
DE.R.crco$gene <- rownames(DE.R.crco)

DE.F.oc2sf.keep <- intersect(DE.F.cgco$gene, DE.F.cgcr$gene) # intersect the DE.F.cgco$gene & DE.F.cgcr$gene
DE.F.oc2sf.exclude <- !(DE.F.oc2sf.keep %in% DE.F.crco$gene) # exclude DE genes in CrCo.
DE.F.MateGene <- DE.F.oc2sf.keep[DE.F.oc2sf.exclude] # perform the exclusion

DE.L.oc2sf.keep <- intersect(DE.L.cgco$gene, DE.L.cgcr$gene)
DE.L.oc2sf.exclude <- !(DE.L.oc2sf.keep %in% DE.L.crco$gene) 
DE.L.MateGene <- DE.L.oc2sf.keep[DE.L.oc2sf.exclude] 

DE.R.oc2sf.keep <- intersect(DE.R.cgco$gene, DE.R.cgcr$gene)
DE.R.oc2sf.exclude <- !(DE.R.oc2sf.keep %in% DE.R.crco$gene) 
DE.R.MateGene <- DE.R.oc2sf.keep[DE.R.oc2sf.exclude] 

message("In flower, number of genes differently expressed between outcrosser and selfer is: ", length(DE.F.MateGene))
message("In leaf, number of genes differently expressed between outcrosser and selfer is: ", length(DE.L.MateGene))
message("In root, number of genes differently expressed between outcrosser and selfer is: ", length(DE.R.MateGene))

# Define genes related with transition from outcrosser to selfer uniqe in flower
# Define 1. Unique DE genes only correlated with mating type transition in flower but not correlated in leaf and root
DE.FlowerUniq.MateGene <- DE.F.MateGene[!(DE.F.MateGene %in% DE.L.MateGene) & !(DE.F.MateGene %in% DE.R.MateGene) ]

# Define 2. Unique MTT related DE genes only expressed in flower but not in others.
# This means all genes have to fit two requirements: 
# 1. Unique DE expressed in Flower
# 2. Related with mating type transition (MTT)

DE.Flower_Union <- union(union(DE.F.cgco$gene, DE.F.cgcr$gene), DE.F.crco$gene) # union genes in Flower
DE.Leaf_Union <- union(union(DE.L.cgco$gene, DE.L.cgcr$gene), DE.L.crco$gene)   # union genes in Leaf
DE.Root_Union <- union(union(DE.R.cgco$gene, DE.R.cgcr$gene), DE.R.crco$gene)   # union genes in Root

DE.Flower.MTT_Uniq <- DE.F.MateGene[!(DE.F.MateGene %in% DE.Leaf_Union) & !(DE.F.MateGene %in% DE.Root_Union) ]
DE.Flower.MTT_Uniq <- cbind(DE.Flower.MTT_Uniq)
DE.Flower.MTT_Uniq <- as.data.frame(DE.Flower.MTT_Uniq)
names(DE.Flower.MTT_Uniq) <- c("Gene.name")

DE.Leaf.MTT_Uniq <- DE.L.MateGene[!(DE.L.MateGene %in% DE.Flower_Union) & !(DE.L.MateGene %in% DE.Root_Union) ]
DE.Leaf.MTT_Uniq <- cbind(DE.Leaf.MTT_Uniq)
DE.Leaf.MTT_Uniq <- as.data.frame(DE.Leaf.MTT_Uniq)
names(DE.Leaf.MTT_Uniq) <- c("Gene.name")
head(DE.Leaf.MTT_Uniq)
dim(DE.Leaf.MTT_Uniq)

DE.Root.MTT_Uniq <- DE.R.MateGene[!(DE.R.MateGene %in% DE.Flower_Union) & !(DE.R.MateGene %in% DE.Leaf_Union) ]
DE.Root.MTT_Uniq <- cbind(DE.Root.MTT_Uniq)
DE.Root.MTT_Uniq <- as.data.frame(DE.Root.MTT_Uniq)
names(DE.Root.MTT_Uniq) <- c("Gene.name")
head(DE.Root.MTT_Uniq)
dim(DE.Root.MTT_Uniq)
# Convert gene names of Capsella to Arabidopsis thaliana 拟南芥
cr2at <- read.table("InputData/AT_CR_mart_export.txt", header=TRUE, sep = "\t")

cr2at <- cr2at[,c("Ortholog.gene.name", "Ortholog.organism.name", "Gene.name", "Organism.name", "Relationship")]
names(cr2at) <- c("Gene.name", "Organism.name", "Ortholog.gene.name", "Ortholog.organism.name", "Relationship")
head(cr2at)
dim(cr2at)

DE.Flower.MTT_Uniq_cr2at <- merge(DE.Flower.MTT_Uniq, cr2at, by = "Gene.name")
dim(DE.Flower.MTT_Uniq_cr2at)

DE.Leaf.MTT_Uniq_cr2at <- merge(DE.Leaf.MTT_Uniq, cr2at, by = "Gene.name")
dim(DE.Leaf.MTT_Uniq_cr2at)

DE.Root.MTT_Uniq_cr2at <- merge(DE.Root.MTT_Uniq, cr2at, by = "Gene.name")
dim(DE.Root.MTT_Uniq_cr2at)

# write table out
# write.table(cbind(DE.Flower.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Flower.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)
# write.table(DE.Flower.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Flower_cr2at.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)
# write.table(cbind(DE.Leaf.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Leaf.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)
# write.table(DE.Leaf.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Leaf_cr2at.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)
# write.table(cbind(DE.Root.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Root.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)
# write.table(DE.Root.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Root_cr2at.txt", append = F, quote = F, 
#             sep = "\t", row.names = F)


length(DE.FlowerUniq.MateGene)
dim(DE.Flower.MTT_Uniq)




### Venn Diagram of pairwise comparison

library(tidyverse) # for the function "alpha" in follow code
#library(hrbrthemes)
#library(tm)
#library(proustr)
#install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(
  x = list(rownames(DE.F.cgco), rownames(DE.F.cgcr), rownames(DE.F.crco)),
  category.names = c("C.g vs C.o" , "C.g vs C.r" , "C.r vs C.o"),
  filename = "EdgeR Pairwise comparison in Flower.png",
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

venn.diagram(
  x = list(rownames(DE.L.cgco), rownames(DE.L.cgcr), rownames(DE.L.crco)),
  category.names = c("C.g vs C.o" , "C.g vs C.r" , "C.r vs C.o"),
  filename = "EdgeR Pairwise comparison in Leaf.png",
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

venn.diagram(
  x = list(rownames(DE.R.cgco), rownames(DE.R.cgcr), rownames(DE.R.crco)),
  category.names = c("C.g vs C.o" , "C.g vs C.r" , "C.r vs C.o"),
  filename = "EdgeR Pairwise comparison in Root.png",
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

venn.diagram(
  x = list(DE.F.MateGene, DE.L.MateGene, DE.R.MateGene),
  category.names = c("Flower" , "Leaf" , "Root"),
  filename = "Flower unique mating type transition relation genes venn plot.png",
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = list(DE.F.MateGene, DE.Leaf_Union, DE.Root_Union),
  category.names = c("MTT genes in Flower" , "Union genes in Leaf" , "Union genes in root"),
  filename = "MTT related DE genes uniqe in flower.png",
  output = TRUE ,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  compression = "lzw",
  lwd = 1,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

















