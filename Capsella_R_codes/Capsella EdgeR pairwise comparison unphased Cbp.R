#############################################################################

#  Pairwise gene expression comparison across Cg, Co, Cr, and unphased Cbp  #

#############################################################################

setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project")

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

#### This file contains reads counts information of each gene in all samples
raw.RC <- read.table("InputData/Cbp_unphased_diploids_counts_FRL_masked_RNA_ASE_unbiased.csv", 
                     header=TRUE, row.names="gene")
str(raw.RC)
head(raw.RC)
dim(raw.RC)
# CPMnorm <- read.table("Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv", header = T, sep = "\t")
# head(CPMnorm)
# dim(CPMnorm)
##### This file contains information of population and tissue in each sample
Accession <- read.csv("InputData/28genomes_unphased_annot_R.csv", stringsAsFactors = F)
head(Accession)
str(Accession)
Accession$species <- factor(Accession$species, levels = c("CG","CR","CO","Cbp"))
### Define populations
######################################################
CR_F <- c(Accession$population=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$population=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$population=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$population=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$population=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$population=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$population=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$population=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$population=="CO" & Accession$tissue=="root")


Cbp_F <- c(Accession$species == "Cbp" & Accession$tissue=="flower")
Cbp_L <- c(Accession$species == "Cbp" & Accession$tissue=="leaf")
Cbp_R <- c(Accession$species == "Cbp" & Accession$tissue=="root")

ASI_F <- c(Accession$population=="ASI" & Accession$tissue=="flower")
EUR_F <- c(Accession$population=="EUR" & Accession$tissue=="flower")
ME_F <- c(Accession$population=="ME" & Accession$tissue=="flower")
CASI_F <- c(Accession$population=="CASI" & Accession$tissue=="flower")

ASI_L <- c(Accession$population=="ASI" & Accession$tissue=="leaf")
EUR_L <- c(Accession$population=="EUR" & Accession$tissue=="leaf")
ME_L <- c(Accession$population=="ME" & Accession$tissue=="leaf")
CASI_L <- c(Accession$population=="CASI" & Accession$tissue=="leaf")

ASI_R <- c(Accession$population=="ASI" & Accession$tissue=="root")
EUR_R <- c(Accession$population=="EUR" & Accession$tissue=="root")
ME_R <- c(Accession$population=="ME" & Accession$tissue=="root")
CASI_R <- c(Accession$population=="CASI" & Accession$tissue=="root")

######################################################
## Flower
F.cgco <- Accession[((Accession$tissue == "flower") & ((Accession$species == "CG") | (Accession$species == "CO"))),]
row.names(F.cgco) <- F.cgco$accession
str(F.cgco)
F.cgcr <- rbind(Accession[(Accession$tissue == "flower") & (Accession$species == "CG"),],
                Accession[(Accession$tissue == "flower") & (Accession$species == "CR"),])
row.names(F.cgcr) <- F.cgcr$accession
F.cgcr
str(F.cgcr)
F.crco <- Accession[((Accession$tissue == "flower") & ((Accession$species == "CR") | (Accession$species == "CO"))),]
row.names(F.crco) <- F.crco$accession
F.crco
str(F.crco)

F.cgcbp <- Accession[((Accession$tissue == "flower") & ((Accession$species == "CG") | (Accession$species == "Cbp"))),]
row.names(F.cgcbp) <- F.cgcbp$accession
str(F.cgcbp)

F.crcbp <- Accession[((Accession$tissue == "flower") & ((Accession$species == "CR") | (Accession$species == "Cbp"))),]
row.names(F.crcbp) <- F.crcbp$accession
str(F.crcbp)

F.cocbp <- Accession[((Accession$tissue == "flower") & ((Accession$species == "CO") | (Accession$species == "Cbp"))),]
row.names(F.cocbp) <- F.cocbp$accession
str(F.cocbp)

## Leaf
L.cgco <- Accession[((Accession$tissue == "leaf") & ((Accession$species == "CG") | (Accession$species == "CO"))),]
row.names(L.cgco) <- L.cgco$accession
str(L.cgco)
L.cgcr <- rbind(Accession[(Accession$tissue == "leaf") & (Accession$species == "CG"),],
                Accession[(Accession$tissue == "leaf") & (Accession$species == "CR"),])
row.names(L.cgcr) <- L.cgcr$accession
str(L.cgcr)
L.crco <- Accession[((Accession$tissue == "leaf") & ((Accession$species == "CR") | (Accession$species == "CO"))),]
row.names(L.crco) <- L.crco$accession
str(L.crco)

L.cgcbp <- Accession[((Accession$tissue == "leaf") & ((Accession$species == "CG") | (Accession$species == "Cbp"))),]
row.names(L.cgcbp) <- L.cgcbp$accession
str(L.cgcbp)

L.crcbp <- Accession[((Accession$tissue == "leaf") & ((Accession$species == "CR") | (Accession$species == "Cbp"))),]
row.names(L.crcbp) <- L.crcbp$accession
str(L.crcbp)

L.cocbp <- Accession[((Accession$tissue == "leaf") & ((Accession$species == "CO") | (Accession$species == "Cbp"))),]
row.names(L.cocbp) <- L.cocbp$accession
str(L.cocbp)

## Root
R.cgco <- Accession[((Accession$tissue == "root") & ((Accession$species == "CG") | (Accession$species == "CO"))),]
row.names(R.cgco) <- R.cgco$accession
str(R.cgco)
R.cgcr <- rbind(Accession[(Accession$tissue == "root") & (Accession$species == "CG"),],
                Accession[(Accession$tissue == "root") & (Accession$species == "CR"),])
row.names(R.cgcr) <- R.cgcr$accession
str(R.cgcr)
R.crco <- Accession[((Accession$tissue == "root") & ((Accession$species == "CR") | (Accession$species == "CO"))),]
row.names(R.crco) <- R.crco$accession
str(R.crco)

R.cgcbp <- Accession[((Accession$tissue == "root") & ((Accession$species == "CG") | (Accession$species == "Cbp"))),]
row.names(R.cgcbp) <- R.cgcbp$accession
str(R.cgcbp)

R.crcbp <- Accession[((Accession$tissue == "root") & ((Accession$species == "CR") | (Accession$species == "Cbp"))),]
row.names(R.crcbp) <- R.crcbp$accession
str(R.crcbp)

R.cocbp <- Accession[((Accession$tissue == "root") & ((Accession$species == "CO") | (Accession$species == "Cbp"))),]
row.names(R.cocbp) <- R.cocbp$accession
str(R.cocbp)
# Filter reads
geneMis <- 1
geneKeep <- cbind(  ((rowSums(raw.RC[,CG_F] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CR_F] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CO_F] >= 10, na.rm = T)) >= geneMis) &
                      
                      ((rowSums(raw.RC[,CG_L] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CR_L] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CO_L] >= 10, na.rm = T)) >= geneMis) &
                      
                      ((rowSums(raw.RC[,CG_R] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CR_R] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CO_R] >= 10, na.rm = T)) >= geneMis) &
                      
                      ((rowSums(raw.RC[,ASI_F] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,EUR_F] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,ME_F] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CASI_F] >= 10, na.rm = T)) >= geneMis) &
                      
                      ((rowSums(raw.RC[,ASI_L] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,EUR_L] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,ME_L] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CASI_L] >= 10, na.rm = T)) >= geneMis) &
                      
                      ((rowSums(raw.RC[,ASI_R] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,EUR_R] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,ME_R] >= 10, na.rm = T)) >= geneMis) &
                      ((rowSums(raw.RC[,CASI_R] >= 10, na.rm = T)) >= geneMis)
)

raw.RC.Keep <- raw.RC[geneKeep,]
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.Keep)[1]))


### DE with EdgeR
# Create EdgeR dataset
# dge = differentical gene expression
F.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_F | CO_F)],         group = F.cgco$species)
F.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_F], raw.RC.Keep[,CR_F]), group = F.cgcr$species)
F.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_F | CO_F)],         group = F.crco$species)
F.cgcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_F | Cbp_F)],         group = F.cgcbp$species)
F.crcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_F | Cbp_F)],         group = F.crcbp$species)
F.cocbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CO_F | Cbp_F)],         group = F.cocbp$species)

L.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_L | CO_L)],         group = L.cgco$species)
L.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_L], raw.RC.Keep[,CR_L]), group = L.cgcr$species)
L.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_L | CO_L)],         group = L.crco$species)
L.cgcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_L | Cbp_L)],         group = L.cgcbp$species)
L.crcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_L | Cbp_L)],         group = L.crcbp$species)
L.cocbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CO_L | Cbp_L)],         group = L.cocbp$species)

R.cgco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_R | CO_R)],         group = R.cgco$species)
R.cgcr.dge <- DGEList(counts = cbind(raw.RC.Keep[,CG_R], raw.RC.Keep[,CR_R]), group = R.cgcr$species)
R.crco.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_R | CO_R)],         group = R.crco$species)
R.cgcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CG_R | Cbp_R)],         group = R.cgcbp$species)
R.crcbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CR_R | Cbp_R)],         group = R.crcbp$species)
R.cocbp.dge <- DGEList(counts = raw.RC.Keep[,rbind(CO_R | Cbp_R)],         group = R.cocbp$species)

# Normalize by total count 
# TMM normalization
F.cgco.dge <- calcNormFactors(F.cgco.dge, method = "TMM")
F.cgcr.dge <- calcNormFactors(F.cgcr.dge, method = "TMM")
F.crco.dge <- calcNormFactors(F.crco.dge, method = "TMM")
F.cgcbp.dge <- calcNormFactors(F.cgcbp.dge, method = "TMM")
F.crcbp.dge <- calcNormFactors(F.crcbp.dge, method = "TMM")
F.cocbp.dge <- calcNormFactors(F.cocbp.dge, method = "TMM")

L.cgco.dge <- calcNormFactors(L.cgco.dge, method = "TMM")
L.cgcr.dge <- calcNormFactors(L.cgcr.dge, method = "TMM")
L.crco.dge <- calcNormFactors(L.crco.dge, method = "TMM")
L.cgcbp.dge <- calcNormFactors(L.cgcbp.dge, method = "TMM")
L.crcbp.dge <- calcNormFactors(L.crcbp.dge, method = "TMM")
L.cocbp.dge <- calcNormFactors(L.cocbp.dge, method = "TMM")

R.cgco.dge <- calcNormFactors(R.cgco.dge, method = "TMM")
R.cgcr.dge <- calcNormFactors(R.cgcr.dge, method = "TMM")
R.crco.dge <- calcNormFactors(R.crco.dge, method = "TMM")
R.cgcbp.dge <- calcNormFactors(R.cgcbp.dge, method = "TMM")
R.crcbp.dge <- calcNormFactors(R.crcbp.dge, method = "TMM")
R.cocbp.dge <- calcNormFactors(R.cocbp.dge, method = "TMM")
# Generate the barplot of library size
barplot(F.cgcbp.dge$samples$lib.size, names.arg=F.cgcbp$accession,
        las =2, border = F, ylab = "Library size (million reads) in Flower", cex.names=0.5)
# Create the contrast matrix
F.cgco.dge.mat <- model.matrix(~ 0 + F.cgco.dge$samples$group)
colnames(F.cgco.dge.mat) <- levels(F.cgco.dge$samples$group)
F.cgcr.dge.mat <- model.matrix(~ 0 + F.cgcr.dge$samples$group)
colnames(F.cgcr.dge.mat) <- levels(F.cgcr.dge$samples$group)
F.crco.dge.mat <- model.matrix(~ 0 + F.crco.dge$samples$group)
colnames(F.crco.dge.mat) <- levels(F.crco.dge$samples$group)
F.cgcbp.dge.mat <- model.matrix(~ 0 + F.cgcbp.dge$samples$group)
colnames(F.cgcbp.dge.mat) <- levels(F.cgcbp.dge$samples$group)
F.crcbp.dge.mat <- model.matrix(~ 0 + F.crcbp.dge$samples$group)
colnames(F.crcbp.dge.mat) <- levels(F.crcbp.dge$samples$group)
F.cocbp.dge.mat <- model.matrix(~ 0 + F.cocbp.dge$samples$group)
colnames(F.cocbp.dge.mat) <- levels(F.cocbp.dge$samples$group)

L.cgco.dge.mat <- model.matrix(~ 0 + L.cgco.dge$samples$group)
colnames(L.cgco.dge.mat) <- levels(L.cgco.dge$samples$group)
L.cgcr.dge.mat <- model.matrix(~ 0 + L.cgcr.dge$samples$group)
colnames(L.cgcr.dge.mat) <- levels(L.cgcr.dge$samples$group)
L.crco.dge.mat <- model.matrix(~ 0 + L.crco.dge$samples$group)
colnames(L.crco.dge.mat) <- levels(L.crco.dge$samples$group)
L.cgcbp.dge.mat <- model.matrix(~ 0 + L.cgcbp.dge$samples$group)
colnames(L.cgcbp.dge.mat) <- levels(L.cgcbp.dge$samples$group)
L.crcbp.dge.mat <- model.matrix(~ 0 + L.crcbp.dge$samples$group)
colnames(L.crcbp.dge.mat) <- levels(L.crcbp.dge$samples$group)
L.cocbp.dge.mat <- model.matrix(~ 0 + L.cocbp.dge$samples$group)
colnames(L.cocbp.dge.mat) <- levels(L.cocbp.dge$samples$group)

R.cgco.dge.mat <- model.matrix(~ 0 + R.cgco.dge$samples$group)
colnames(R.cgco.dge.mat) <- levels(R.cgco.dge$samples$group)
R.cgcr.dge.mat <- model.matrix(~ 0 + R.cgcr.dge$samples$group)
colnames(R.cgcr.dge.mat) <- levels(R.cgcr.dge$samples$group)
R.crco.dge.mat <- model.matrix(~ 0 + R.crco.dge$samples$group)
colnames(R.crco.dge.mat) <- levels(R.crco.dge$samples$group)
R.cgcbp.dge.mat <- model.matrix(~ 0 + R.cgcbp.dge$samples$group)
colnames(R.cgcbp.dge.mat) <- levels(R.cgcbp.dge$samples$group)
R.crcbp.dge.mat <- model.matrix(~ 0 + R.crcbp.dge$samples$group)
colnames(R.crcbp.dge.mat) <- levels(R.crcbp.dge$samples$group)
R.cocbp.dge.mat <- model.matrix(~ 0 + R.cocbp.dge$samples$group)
colnames(R.cocbp.dge.mat) <- levels(R.cocbp.dge$samples$group)

# Estimate dispersion parameter for GLM
# Flower
F.cgco.dge <- estimateGLMCommonDisp(F.cgco.dge, F.cgco.dge.mat)
F.cgco.dge <- estimateGLMTrendedDisp(F.cgco.dge, F.cgco.dge.mat, method="power")
F.cgco.dge <- estimateGLMTagwiseDisp(F.cgco.dge, F.cgco.dge.mat)

F.cgcr.dge <- estimateGLMCommonDisp(F.cgcr.dge, F.cgcr.dge.mat)
F.cgcr.dge <- estimateGLMTrendedDisp(F.cgcr.dge, F.cgcr.dge.mat, method="power")
F.cgcr.dge <- estimateGLMTagwiseDisp(F.cgcr.dge, F.cgcr.dge.mat)

F.crco.dge <- estimateGLMCommonDisp(F.crco.dge, F.crco.dge.mat)
F.crco.dge <- estimateGLMTrendedDisp(F.crco.dge, F.crco.dge.mat, method="power")
F.crco.dge <- estimateGLMTagwiseDisp(F.crco.dge, F.crco.dge.mat)

F.cgcbp.dge <- estimateGLMCommonDisp(F.cgcbp.dge, F.cgcbp.dge.mat)
F.cgcbp.dge <- estimateGLMTrendedDisp(F.cgcbp.dge, F.cgcbp.dge.mat, method="power")
F.cgcbp.dge <- estimateGLMTagwiseDisp(F.cgcbp.dge, F.cgcbp.dge.mat)

F.crcbp.dge <- estimateGLMCommonDisp(F.crcbp.dge, F.crcbp.dge.mat)
F.crcbp.dge <- estimateGLMTrendedDisp(F.crcbp.dge, F.crcbp.dge.mat, method="power")
F.crcbp.dge <- estimateGLMTagwiseDisp(F.crcbp.dge, F.crcbp.dge.mat)

F.cocbp.dge <- estimateGLMCommonDisp(F.cocbp.dge, F.cocbp.dge.mat)
F.cocbp.dge <- estimateGLMTrendedDisp(F.cocbp.dge, F.cocbp.dge.mat, method="power")
F.cocbp.dge <- estimateGLMTagwiseDisp(F.cocbp.dge, F.cocbp.dge.mat)

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

L.cgcbp.dge <- estimateGLMCommonDisp(L.cgcbp.dge, L.cgcbp.dge.mat)
L.cgcbp.dge <- estimateGLMTrendedDisp(L.cgcbp.dge, L.cgcbp.dge.mat, method="power")
L.cgcbp.dge <- estimateGLMTagwiseDisp(L.cgcbp.dge, L.cgcbp.dge.mat)

L.crcbp.dge <- estimateGLMCommonDisp(L.crcbp.dge, L.crcbp.dge.mat)
L.crcbp.dge <- estimateGLMTrendedDisp(L.crcbp.dge, L.crcbp.dge.mat, method="power")
L.crcbp.dge <- estimateGLMTagwiseDisp(L.crcbp.dge, L.crcbp.dge.mat)

L.cocbp.dge <- estimateGLMCommonDisp(L.cocbp.dge, L.cocbp.dge.mat)
L.cocbp.dge <- estimateGLMTrendedDisp(L.cocbp.dge, L.cocbp.dge.mat, method="power")
L.cocbp.dge <- estimateGLMTagwiseDisp(L.cocbp.dge, L.cocbp.dge.mat)

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

R.cgcbp.dge <- estimateGLMCommonDisp(R.cgcbp.dge, R.cgcbp.dge.mat)
R.cgcbp.dge <- estimateGLMTrendedDisp(R.cgcbp.dge, R.cgcbp.dge.mat, method="power")
R.cgcbp.dge <- estimateGLMTagwiseDisp(R.cgcbp.dge, R.cgcbp.dge.mat)

R.crcbp.dge <- estimateGLMCommonDisp(R.crcbp.dge, R.crcbp.dge.mat)
R.crcbp.dge <- estimateGLMTrendedDisp(R.crcbp.dge, R.crcbp.dge.mat, method="power")
R.crcbp.dge <- estimateGLMTagwiseDisp(R.crcbp.dge, R.crcbp.dge.mat)

R.cocbp.dge <- estimateGLMCommonDisp(R.cocbp.dge, R.cocbp.dge.mat)
R.cocbp.dge <- estimateGLMTrendedDisp(R.cocbp.dge, R.cocbp.dge.mat, method="power")
R.cocbp.dge <- estimateGLMTagwiseDisp(R.cocbp.dge, R.cocbp.dge.mat)
# Plot mean-variance
plotBCV(F.cgco.dge)
plotBCV(F.cgcr.dge)
plotBCV(F.crco.dge)
plotBCV(F.cgcbp.dge)
plotBCV(F.crcbp.dge)
plotBCV(F.cocbp.dge)
plotBCV(L.cgco.dge)
plotBCV(L.cgcr.dge)
plotBCV(L.crco.dge)
plotBCV(L.cgcbp.dge)
plotBCV(L.crcbp.dge)
plotBCV(L.cocbp.dge)
plotBCV(R.cgco.dge)
plotBCV(R.cgcr.dge)
plotBCV(R.crco.dge)
plotBCV(R.cgcbp.dge)
plotBCV(R.crcbp.dge)
plotBCV(R.cocbp.dge)
##########
## DEG analysis 
# Set threshold
p.threshold <- 0.01
## edgeR ##
# Model fitting
F.cgco.fit <- glmFit(F.cgco.dge, F.cgco.dge.mat)
F.cgcr.fit <- glmFit(F.cgcr.dge, F.cgcr.dge.mat)
F.crco.fit <- glmFit(F.crco.dge, F.crco.dge.mat)
F.cgcbp.fit <- glmFit(F.cgcbp.dge, F.cgcbp.dge.mat)
F.crcbp.fit <- glmFit(F.crcbp.dge, F.crcbp.dge.mat)
F.cocbp.fit <- glmFit(F.cocbp.dge, F.cocbp.dge.mat)

L.cgco.fit <- glmFit(L.cgco.dge, L.cgco.dge.mat)
L.cgcr.fit <- glmFit(L.cgcr.dge, L.cgcr.dge.mat)
L.crco.fit <- glmFit(L.crco.dge, L.crco.dge.mat)
L.cgcbp.fit <- glmFit(L.cgcbp.dge, L.cgcbp.dge.mat)
L.crcbp.fit <- glmFit(L.crcbp.dge, L.crcbp.dge.mat)
L.cocbp.fit <- glmFit(L.cocbp.dge, L.cocbp.dge.mat)

R.cgco.fit <- glmFit(R.cgco.dge, R.cgco.dge.mat)
R.cgcr.fit <- glmFit(R.cgcr.dge, R.cgcr.dge.mat)
R.crco.fit <- glmFit(R.crco.dge, R.crco.dge.mat)
R.cgcbp.fit <- glmFit(R.cgcbp.dge, R.cgcbp.dge.mat)
R.crcbp.fit <- glmFit(R.crcbp.dge, R.crcbp.dge.mat)
R.cocbp.fit <- glmFit(R.cocbp.dge, R.cocbp.dge.mat)

# Differential expression
# Flower
F.cgco.contrast <- makeContrasts(CG - CO, levels = F.cgco.dge.mat)
F.cgco.lrt <- glmLRT(F.cgco.fit, contrast = F.cgco.contrast)

F.cgcr.contrast <- makeContrasts(CG - CR, levels = F.cgcr.dge.mat)
F.cgcr.lrt <- glmLRT(F.cgcr.fit, contrast = F.cgcr.contrast)

F.crco.contrast <- makeContrasts(CR - CO, levels = F.crco.dge.mat)
F.crco.lrt <- glmLRT(F.crco.fit, contrast = F.crco.contrast)

F.cgcbp.contrast <- makeContrasts(CG - Cbp, levels = F.cgcbp.dge.mat)
F.cgcbp.lrt <- glmLRT(F.cgcbp.fit, contrast = F.cgcbp.contrast)

F.crcbp.contrast <- makeContrasts(CR - Cbp, levels = F.crcbp.dge.mat)
F.crcbp.lrt <- glmLRT(F.crcbp.fit, contrast = F.crcbp.contrast)

F.cocbp.contrast <- makeContrasts(CO - Cbp, levels = F.cocbp.dge.mat)
F.cocbp.lrt <- glmLRT(F.cocbp.fit, contrast = F.cocbp.contrast)
# Leaf
L.cgco.contrast <- makeContrasts(CG - CO, levels = L.cgco.dge.mat)
L.cgco.lrt <- glmLRT(L.cgco.fit, contrast = L.cgco.contrast)

L.cgcr.contrast <- makeContrasts(CG - CR, levels = L.cgcr.dge.mat)
L.cgcr.lrt <- glmLRT(L.cgcr.fit, contrast = L.cgcr.contrast)

L.crco.contrast <- makeContrasts(CR - CO, levels = L.crco.dge.mat)
L.crco.lrt <- glmLRT(L.crco.fit, contrast = L.crco.contrast)

L.cgcbp.contrast <- makeContrasts(CG - Cbp, levels = L.cgcbp.dge.mat)
L.cgcbp.lrt <- glmLRT(L.cgcbp.fit, contrast = L.cgcbp.contrast)

L.crcbp.contrast <- makeContrasts(CR - Cbp, levels = L.crcbp.dge.mat)
L.crcbp.lrt <- glmLRT(L.crcbp.fit, contrast = L.crcbp.contrast)

L.cocbp.contrast <- makeContrasts(CO - Cbp, levels = L.cocbp.dge.mat)
L.cocbp.lrt <- glmLRT(L.cocbp.fit, contrast = L.cocbp.contrast)
# Root
R.cgco.contrast <- makeContrasts(CG - CO, levels = R.cgco.dge.mat)
R.cgco.lrt <- glmLRT(R.cgco.fit, contrast = R.cgco.contrast)

R.cgcr.contrast <- makeContrasts(CG - CR, levels = R.cgcr.dge.mat)
R.cgcr.lrt <- glmLRT(R.cgcr.fit, contrast = R.cgcr.contrast)

R.crco.contrast <- makeContrasts(CR - CO, levels = R.crco.dge.mat)
R.crco.lrt <- glmLRT(R.crco.fit, contrast = R.crco.contrast)

R.cgcbp.contrast <- makeContrasts(CG - Cbp, levels = R.cgcbp.dge.mat)
R.cgcbp.lrt <- glmLRT(R.cgcbp.fit, contrast = R.cgcbp.contrast)

R.crcbp.contrast <- makeContrasts(CR - Cbp, levels = R.crcbp.dge.mat)
R.crcbp.lrt <- glmLRT(R.crcbp.fit, contrast = R.crcbp.contrast)

R.cocbp.contrast <- makeContrasts(CO - Cbp, levels = R.cocbp.dge.mat)
R.cocbp.lrt <- glmLRT(R.cocbp.fit, contrast = R.cocbp.contrast)

# Access results tables
F.cgco.table <- F.cgco.lrt$table
F.cgcr.table <- F.cgcr.lrt$table
F.crco.table <- F.crco.lrt$table
F.cgcbp.table <- F.cgcbp.lrt$table
F.crcbp.table <- F.crcbp.lrt$table
F.cocbp.table <- F.cocbp.lrt$table

L.cgco.table <- L.cgco.lrt$table
L.cgcr.table <- L.cgcr.lrt$table
L.crco.table <- L.crco.lrt$table
L.cgcbp.table <- L.cgcbp.lrt$table
L.crcbp.table <- L.crcbp.lrt$table
L.cocbp.table <- L.cocbp.lrt$table

R.cgco.table <- R.cgco.lrt$table
R.cgcr.table <- R.cgcr.lrt$table
R.crco.table <- R.crco.lrt$table
R.cgcbp.table <- R.cgcbp.lrt$table
R.crcbp.table <- R.crcbp.lrt$table
R.cocbp.table <- R.cocbp.lrt$table

# Toptags Result with FDR
F.cgco.topTags <- topTags(F.cgco.lrt, n = nrow(F.cgco.lrt$table))$table
F.cgcr.topTags <- topTags(F.cgcr.lrt, n = nrow(F.cgcr.lrt$table))$table
F.crco.topTags <- topTags(F.crco.lrt, n = nrow(F.crco.lrt$table))$table
F.cgcbp.topTags <- topTags(F.cgcbp.lrt, n = nrow(F.cgcbp.lrt$table))$table
F.crcbp.topTags <- topTags(F.crcbp.lrt, n = nrow(F.crcbp.lrt$table))$table
F.cocbp.topTags <- topTags(F.cocbp.lrt, n = nrow(F.cocbp.lrt$table))$table

L.cgco.topTags <- topTags(L.cgco.lrt, n = nrow(L.cgco.lrt$table))$table
L.cgcr.topTags <- topTags(L.cgcr.lrt, n = nrow(L.cgcr.lrt$table))$table
L.crco.topTags <- topTags(L.crco.lrt, n = nrow(L.crco.lrt$table))$table
L.cgcbp.topTags <- topTags(L.cgcbp.lrt, n = nrow(L.cgcbp.lrt$table))$table
L.crcbp.topTags <- topTags(L.crcbp.lrt, n = nrow(L.crcbp.lrt$table))$table
L.cocbp.topTags <- topTags(L.cocbp.lrt, n = nrow(L.cocbp.lrt$table))$table

R.cgco.topTags <- topTags(R.cgco.lrt, n = nrow(R.cgco.lrt$table))$table
R.cgcr.topTags <- topTags(R.cgcr.lrt, n = nrow(R.cgcr.lrt$table))$table
R.crco.topTags <- topTags(R.crco.lrt, n = nrow(R.crco.lrt$table))$table
R.cgcbp.topTags <- topTags(R.cgcbp.lrt, n = nrow(R.cgcbp.lrt$table))$table
R.crcbp.topTags <- topTags(R.crcbp.lrt, n = nrow(R.crcbp.lrt$table))$table
R.cocbp.topTags <- topTags(R.cocbp.lrt, n = nrow(R.cocbp.lrt$table))$table

# Extract the DE genes with FDR < 0.01
DE.F.cgco <- F.cgco.topTags[(abs(F.cgco.topTags$logFC) > 1) & (F.cgco.topTags$FDR <= 0.01),]
DE.F.cgcr <- F.cgcr.topTags[(abs(F.cgcr.topTags$logFC) > 1) & (F.cgcr.topTags$FDR <= 0.01),]
DE.F.crco <- F.crco.topTags[(abs(F.crco.topTags$logFC) > 1) & (F.crco.topTags$FDR <= 0.01),]
DE.F.cgcbp <- F.cgcbp.topTags[(abs(F.cgcbp.topTags$logFC) > 1) & (F.cgcbp.topTags$FDR <= 0.01),]
DE.F.crcbp <- F.crcbp.topTags[(abs(F.crcbp.topTags$logFC) > 1) & (F.crcbp.topTags$FDR <= 0.01),]
DE.F.cocbp <- F.cocbp.topTags[(abs(F.cocbp.topTags$logFC) > 1) & (F.cocbp.topTags$FDR <= 0.01),]

DE.L.cgco <- L.cgco.topTags[(abs(L.cgco.topTags$logFC) > 1) & (L.cgco.topTags$FDR <= 0.01),]
DE.L.cgcr <- L.cgcr.topTags[(abs(L.cgcr.topTags$logFC) > 1) & (L.cgcr.topTags$FDR <= 0.01),]
DE.L.crco <- L.crco.topTags[(abs(L.crco.topTags$logFC) > 1) & (L.crco.topTags$FDR <= 0.01),]
DE.L.cgcbp <- L.cgcbp.topTags[(abs(L.cgcbp.topTags$logFC) > 1) & (L.cgcbp.topTags$FDR <= 0.01),]
DE.L.crcbp <- L.crcbp.topTags[(abs(L.crcbp.topTags$logFC) > 1) & (L.crcbp.topTags$FDR <= 0.01),]
DE.L.cocbp <- L.cocbp.topTags[(abs(L.cocbp.topTags$logFC) > 1) & (L.cocbp.topTags$FDR <= 0.01),]

DE.R.cgco <- R.cgco.topTags[(abs(R.cgco.topTags$logFC) > 1) & (R.cgco.topTags$FDR <= 0.01),]
DE.R.cgcr <- R.cgcr.topTags[(abs(R.cgcr.topTags$logFC) > 1) & (R.cgcr.topTags$FDR <= 0.01),]
DE.R.crco <- R.crco.topTags[(abs(R.crco.topTags$logFC) > 1) & (R.crco.topTags$FDR <= 0.01),]
DE.R.cgcbp <- R.cgcbp.topTags[(abs(R.cgcbp.topTags$logFC) > 1) & (R.cgcbp.topTags$FDR <= 0.01),]
DE.R.crcbp <- R.crcbp.topTags[(abs(R.crcbp.topTags$logFC) > 1) & (R.crcbp.topTags$FDR <= 0.01),]
DE.R.cocbp <- R.cocbp.topTags[(abs(R.cocbp.topTags$logFC) > 1) & (R.cocbp.topTags$FDR <= 0.01),]

message("Number of DE genes between Cg and Co in Flower ", dim(DE.F.cgco)[1])
message("Number of DE genes between Cg and Cr in Flower ", dim(DE.F.cgcr)[1])
message("Number of DE genes between Cr and Co in Flower ", dim(DE.F.crco)[1])
message("Number of DE genes between Cg and Cbp in Flower ", dim(DE.F.cgcbp)[1])
message("Number of DE genes between Cr and Cbp in Flower ", dim(DE.F.crcbp)[1])
message("Number of DE genes between Co and Cbp in Flower ", dim(DE.F.cocbp)[1])

message("Number of DE genes between Cg and Co in Leaf ", dim(DE.L.cgco)[1])
message("Number of DE genes between Cg and Cr in Leaf ", dim(DE.L.cgcr)[1])
message("Number of DE genes between Cr and Co in Leaf ", dim(DE.L.crco)[1])
message("Number of DE genes between Cg and Cbp in Leaf ", dim(DE.L.cgcbp)[1])
message("Number of DE genes between Cr and Cbp in Leaf ", dim(DE.L.crcbp)[1])
message("Number of DE genes between Co and Cbp in Leaf ", dim(DE.L.cocbp)[1])

message("Number of DE genes between Cg and Co in Root ", dim(DE.R.cgco)[1])
message("Number of DE genes between Cg and Cr in Root ", dim(DE.R.cgcr)[1])
message("Number of DE genes between Cr and Co in Root ", dim(DE.R.crco)[1])
message("Number of DE genes between Cg and Cbp in Root ", dim(DE.R.cgcbp)[1])
message("Number of DE genes between Cr and Cbp in Root ", dim(DE.R.crcbp)[1])
message("Number of DE genes between Co and Cbp in Root ", dim(DE.R.cocbp)[1])
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
# Need to add Cbp
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
     main="Gene expression changes in Cg Co within flower", xlab="logFC", ylab = "Counts")
hist(F.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cg Cr within flower", xlab="logFC", ylab = "Counts")
hist(F.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cr Co within flower", xlab="logFC", ylab = "Counts")

hist(L.cgco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cg Co within leaf", xlab="logFC", ylab = "Counts")
hist(L.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cg Cr within leaf", xlab="logFC", ylab = "Counts")
hist(L.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cr Co within leaf", xlab="logFC", ylab = "Counts")

hist(R.cgco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cg Co within root", xlab="logFC", ylab = "Counts")
hist(R.cgcr.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cg Cr within root", xlab="logFC", ylab = "Counts")
hist(R.crco.topTags$logFC, breaks = 100, xlim = c(-15,20), ylim = c(0,5000), 
     main="Gene expression changes in Cr Co within root", xlab="logFC", ylab = "Counts")

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
                          ggtheme = theme_bw()) +
  scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
  stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") + # add mean value
  rotate()  + # rotate the boxplot 
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
                          ggtheme = theme_bw()) +
  scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
  stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
  rotate()  + rremove("legend")

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
                          alpha = 0.5, ggtheme = theme_bw()) +
  scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
  stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") +
  rotate()  + rremove("legend")

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
write.table(cbind(DE.Flower.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Flower.txt", append = F, quote = F, 
            sep = "\t", row.names = F)
write.table(DE.Flower.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Flower_cr2at.txt", append = F, quote = F, 
            sep = "\t", row.names = F)
write.table(cbind(DE.Leaf.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Leaf.txt", append = F, quote = F, 
            sep = "\t", row.names = F)
write.table(DE.Leaf.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Leaf_cr2at.txt", append = F, quote = F, 
            sep = "\t", row.names = F)
write.table(cbind(DE.Root.MTT_Uniq), file = "OutputData/DE_genes_of_MTT_unique_in_Root.txt", append = F, quote = F, 
            sep = "\t", row.names = F)
write.table(DE.Root.MTT_Uniq_cr2at, file = "OutputData/DE_genes_of_MTT_unique_in_Root_cr2at.txt", append = F, quote = F, 
            sep = "\t", row.names = F)


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


#############################################################################################
# Levels of gene expression in C. bursa-pastoris relative to its parental species Cr and Co.

# check the data format and dimension
head(F.crcbp.topTags)
dim(F.crcbp.topTags)

# add gene name column to dataframe
F.crcbp.topTags$Gene.name <- row.names(F.crcbp.topTags) # crcbp in flower
F.cocbp.topTags$Gene.name <- row.names(F.cocbp.topTags) # cocbp in flower
F.crco.topTags$Gene.name  <- row.names(F.crco.topTags)  # crco  in flower

L.crcbp.topTags$Gene.name <- row.names(L.crcbp.topTags) # crcbp in leaf
L.cocbp.topTags$Gene.name <- row.names(L.cocbp.topTags) # cocbp in leaf
L.crco.topTags$Gene.name  <- row.names(L.crco.topTags)  # crco  in leaf

R.crcbp.topTags$Gene.name <- row.names(R.crcbp.topTags) # crcbp in root
R.cocbp.topTags$Gene.name <- row.names(R.cocbp.topTags) # cocbp in root
R.crco.topTags$Gene.name  <- row.names(R.crco.topTags)  # crco  in root

# Sort each dataframe by gene name
F_crcbp_sort <- F.crcbp.topTags[order(F.crcbp.topTags$Gene.name),]
F_cocbp_sort <- F.cocbp.topTags[order(F.cocbp.topTags$Gene.name),]
F_crco_sort  <- F.crco.topTags[order(F.crco.topTags$Gene.name),]

L_crcbp_sort <- L.crcbp.topTags[order(L.crcbp.topTags$Gene.name),]
L_cocbp_sort <- L.cocbp.topTags[order(L.cocbp.topTags$Gene.name),]
L_crco_sort  <- L.crco.topTags[order(L.crco.topTags$Gene.name),]

R_crcbp_sort <- R.crcbp.topTags[order(R.crcbp.topTags$Gene.name),]
R_cocbp_sort <- R.cocbp.topTags[order(R.cocbp.topTags$Gene.name),]
R_crco_sort  <- R.crco.topTags[order(R.crco.topTags$Gene.name),]

# check the sorted result
head(F_crcbp_sort) 
head(F_cocbp_sort)
head(F_crco_sort)
head(L_crcbp_sort)
head(L_cocbp_sort)
head(L_crco_sort)
head(R_crcbp_sort)
head(R_cocbp_sort)
head(R_crco_sort)

# combine three comparison to a single file
F_CbpCrCo <- cbind(F_crco_sort[,c("Gene.name","logFC","FDR")], 
                   F_crcbp_sort[,c("logFC","FDR")], 
                   F_cocbp_sort[,c("logFC","FDR")])
names(F_CbpCrCo) <- c("Gene.name","F_crco_logFC","F_crco_FDR",
                      "F_crcbp_logFC","F_crcbp_FDR",
                      "F_cocbp_logFC","F_cocbp_FDR")
head(F_CbpCrCo)

L_CbpCrCo <- cbind(L_crco_sort[,c("Gene.name","logFC","FDR")], 
                   L_crcbp_sort[,c("logFC","FDR")], 
                   L_cocbp_sort[,c("logFC","FDR")])
names(L_CbpCrCo) <- c("Gene.name","L_crco_logFC","L_crco_FDR",
                      "L_crcbp_logFC","L_crcbp_FDR",
                      "L_cocbp_logFC","L_cocbp_FDR")
head(L_CbpCrCo)

R_CbpCrCo <- cbind(R_crco_sort[,c("Gene.name","logFC","FDR")], 
                   R_crcbp_sort[,c("logFC","FDR")], 
                   R_cocbp_sort[,c("logFC","FDR")])
names(R_CbpCrCo) <- c("Gene.name","R_crco_logFC","R_crco_FDR",
                      "R_crcbp_logFC","R_crcbp_FDR",
                      "R_cocbp_logFC","R_cocbp_FDR")
head(R_CbpCrCo)

## Flower
# Expression pattern of No Difference: Cr = Cbp = Co
F_ND <- c((F_CbpCrCo$F_crcbp_FDR >= 0.05) &
            (F_CbpCrCo$F_cocbp_FDR >= 0.05))

F_NoDifference <- F_CbpCrCo[F_ND,]
message("In flower, number of genes with expression pattern of No Difference is: ", dim(F_NoDifference)[1], 
        ", (", round(dim(F_NoDifference)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
F_Ad_1 <- c(# (F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
  ((F_CbpCrCo$F_crcbp_logFC > 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
    ((F_CbpCrCo$F_cocbp_logFC < 0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_Ad_2 <- c(# (F_CbpCrCo$F_crco_logFC < -1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
  ((F_CbpCrCo$F_crcbp_logFC < 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
    ((F_CbpCrCo$F_cocbp_logFC > 0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_Additivity <- F_CbpCrCo[F_Ad_1 | F_Ad_2,]
message("In flower, number of genes with expression pattern of Additivity is: ", dim(F_Additivity)[1],
        ", (", round(dim(F_Additivity)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
F_Dom_1 <- c( #(F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
  (F_CbpCrCo$F_crcbp_FDR >= 0.05) &
    ((F_CbpCrCo$F_cocbp_logFC < -0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_Dom_2 <- c( #F_CbpCrCo$F_crco_logFC < -1 & F_CbpCrCo$F_crco_FDR < 0.05 &
  (F_CbpCrCo$F_crcbp_FDR >= 0.05) &
    ((F_CbpCrCo$F_cocbp_logFC > 0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_CrDom <- F_CbpCrCo[F_Dom_1 | F_Dom_2,]
message("In flower, number of genes with expression pattern of Cr Dominance is: ", dim(F_CrDom)[1],
        ", (", round(dim(F_CrDom)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
F_Dom_3 <- c( #(F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
  ((F_CbpCrCo$F_crcbp_logFC > 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
    (F_CbpCrCo$F_cocbp_FDR  >= 0.05))
F_Dom_4 <- c( #(F_CbpCrCo$F_crco_logFC < -1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
  ((F_CbpCrCo$F_crcbp_logFC < 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
    (F_CbpCrCo$F_cocbp_FDR >= 0.05))
F_CoDom <- F_CbpCrCo[F_Dom_3 | F_Dom_4,]
message("In flower, number of genes with expression pattern of Co Dominance is: ", dim(F_CoDom)[1],
        ", (", round(dim(F_CoDom)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
F_Dominance <- F_CbpCrCo[F_Dom_1 | F_Dom_2 | F_Dom_3 | F_Dom_4,]
message("In flower, number of genes with expression pattern of Dominance is: ", dim(F_Dominance)[1],
        ", (", round(dim(F_Dominance)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
F_PTS <- c(((F_CbpCrCo$F_crcbp_logFC < 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
             ((F_CbpCrCo$F_cocbp_logFC < 0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_PositiveTransgressive <- F_CbpCrCo[F_PTS,]
message("In flower, number of genes with expression pattern of Positive Transgressive is: ", dim(F_PositiveTransgressive)[1],
        ", (", round(dim(F_PositiveTransgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
F_NTS <- c(((F_CbpCrCo$F_crcbp_logFC > 0) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
             ((F_CbpCrCo$F_cocbp_logFC > 0) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
F_NegativeTransgressive <- F_CbpCrCo[F_NTS,]
message("In flower, number of genes with expression pattern of Negative Transgressive is: ", dim(F_NegativeTransgressive)[1],
        ", (", round(dim(F_NegativeTransgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
F_Transgressive <- F_CbpCrCo[F_PTS | F_NTS, ]
message("In flower, number of genes with expression pattern of Negative Transgressive is: ", dim(F_Transgressive)[1],
        ", (", round(dim(F_Transgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

## Leaf
# Expression pattern of No Difference: Cr = Cbp = Co
L_ND <- c((L_CbpCrCo$L_crcbp_FDR >= 0.05) &
            (L_CbpCrCo$L_cocbp_FDR >= 0.05))

L_NoDifference <- L_CbpCrCo[L_ND,]
message("In leaf, number of genes with expression pattern of No Difference is: ", dim(L_NoDifference)[1], 
        ", (", round(dim(L_NoDifference)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
L_Ad_1 <- c(# (L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC > 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC < 0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_Ad_2 <- c(# (L_CbpCrCo$L_crco_logFC < -1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC < 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC > 0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_Additivity <- L_CbpCrCo[L_Ad_1 | L_Ad_2,]
message("In leaf, number of genes with expression pattern of Additivity is: ", dim(L_Additivity)[1],
        ", (", round(dim(L_Additivity)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
L_Dom_1 <- c( #(L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  (L_CbpCrCo$L_crcbp_FDR >= 0.05) &
    ((L_CbpCrCo$L_cocbp_logFC < -0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_Dom_2 <- c( #L_CbpCrCo$L_crco_logFC < -1 & L_CbpCrCo$L_crco_FDR < 0.05 &
  (L_CbpCrCo$L_crcbp_FDR >= 0.05) &
    ((L_CbpCrCo$L_cocbp_logFC > 0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_CrDom <- L_CbpCrCo[L_Dom_1 | L_Dom_2,]
message("In leaf, number of genes with expression pattern of Cr Dominance is: ", dim(L_CrDom)[1],
        ", (", round(dim(L_CrDom)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
L_Dom_3 <- c( #(L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC > 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    (L_CbpCrCo$L_cocbp_FDR  >= 0.05))
L_Dom_4 <- c( #(L_CbpCrCo$L_crco_logFC < -1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC < 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    (L_CbpCrCo$L_cocbp_FDR >= 0.05))
L_CoDom <- L_CbpCrCo[L_Dom_3 | L_Dom_4,]
message("In leaf, number of genes with expression pattern of Co Dominance is: ", dim(L_CoDom)[1],
        ", (", round(dim(L_CoDom)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
L_Dominance <- L_CbpCrCo[L_Dom_1 | L_Dom_2 | L_Dom_3 | L_Dom_4,]
message("In leaf, number of genes with expression pattern of Dominance is: ", dim(L_Dominance)[1],
        ", (", round(dim(L_Dominance)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
L_PTS <- c(((L_CbpCrCo$L_crcbp_logFC < 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
             ((L_CbpCrCo$L_cocbp_logFC < 0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_PositiveTransgressive <- L_CbpCrCo[L_PTS,]
message("In leaf, number of genes with expression pattern of Positive Transgressive is: ", dim(L_PositiveTransgressive)[1],
        ", (", round(dim(L_PositiveTransgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
L_NTS <- c(((L_CbpCrCo$L_crcbp_logFC > 0) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
             ((L_CbpCrCo$L_cocbp_logFC > 0) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
L_NegativeTransgressive <- L_CbpCrCo[L_NTS,]
message("In leaf, number of genes with expression pattern of Negative Transgressive is: ", dim(L_NegativeTransgressive)[1],
        ", (", round(dim(L_NegativeTransgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
L_Transgressive <- L_CbpCrCo[L_PTS | L_NTS, ]
message("In leaf, number of genes with expression pattern of Negative Transgressive is: ", dim(L_Transgressive)[1],
        ", (", round(dim(L_Transgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

## Root
# Expression pattern of No Difference: Cr = Cbp = Co
R_ND <- c((R_CbpCrCo$R_crcbp_FDR >= 0.05) &
            (R_CbpCrCo$R_cocbp_FDR >= 0.05))

R_NoDifference <- R_CbpCrCo[R_ND,]
message("In root, number of genes with expression pattern of No Difference is: ", dim(R_NoDifference)[1], 
        ", (", round(dim(R_NoDifference)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
R_Ad_1 <- c(# (R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC > 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC < 0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_Ad_2 <- c(# (R_CbpCrCo$R_crco_logFC < -1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC < 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC > 0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_Additivity <- R_CbpCrCo[R_Ad_1 | R_Ad_2,]
message("In root, number of genes with expression pattern of Additivity is: ", dim(R_Additivity)[1],
        ", (", round(dim(R_Additivity)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
R_Dom_1 <- c( #(R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  (R_CbpCrCo$R_crcbp_FDR >= 0.05) &
    ((R_CbpCrCo$R_cocbp_logFC < -0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_Dom_2 <- c( #R_CbpCrCo$R_crco_logFC < -1 & R_CbpCrCo$R_crco_FDR < 0.05 &
  (R_CbpCrCo$R_crcbp_FDR >= 0.05) &
    ((R_CbpCrCo$R_cocbp_logFC > 0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_CrDom <- R_CbpCrCo[R_Dom_1 | R_Dom_2,]
message("In root, number of genes with expression pattern of Cr Dominance is: ", dim(R_CrDom)[1],
        ", (", round(dim(R_CrDom)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
R_Dom_3 <- c( #(R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC > 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    (R_CbpCrCo$R_cocbp_FDR  >= 0.05))
R_Dom_4 <- c( #(R_CbpCrCo$R_crco_logFC < -1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC < 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    (R_CbpCrCo$R_cocbp_FDR >= 0.05))
R_CoDom <- R_CbpCrCo[R_Dom_3 | R_Dom_4,]
message("In root, number of genes with expression pattern of Co Dominance is: ", dim(R_CoDom)[1],
        ", (", round(dim(R_CoDom)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
R_Dominance <- R_CbpCrCo[R_Dom_1 | R_Dom_2 | R_Dom_3 | R_Dom_4,]
message("In root, number of genes with expression pattern of Dominance is: ", dim(R_Dominance)[1],
        ", (", round(dim(R_Dominance)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
R_PTS <- c(((R_CbpCrCo$R_crcbp_logFC < 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
             ((R_CbpCrCo$R_cocbp_logFC < 0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_PositiveTransgressive <- R_CbpCrCo[R_PTS,]
message("In root, number of genes with expression pattern of Positive Transgressive is: ", dim(R_PositiveTransgressive)[1],
        ", (", round(dim(R_PositiveTransgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
R_NTS <- c(((R_CbpCrCo$R_crcbp_logFC > 0) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
             ((R_CbpCrCo$R_cocbp_logFC > 0) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
R_NegativeTransgressive <- R_CbpCrCo[R_NTS,]
message("In root, number of genes with expression pattern of Negative Transgressive is: ", dim(R_NegativeTransgressive)[1],
        ", (", round(dim(R_NegativeTransgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
R_Transgressive <- R_CbpCrCo[R_PTS | R_NTS, ]
message("In root, number of genes with expression pattern of Negative Transgressive is: ", dim(R_Transgressive)[1],
        ", (", round(dim(R_Transgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")





################################################################################
## Flower, Strict, Fold chang > 2;
# Expression pattern of No Difference: Cr = Cbp = Co
FS_ND <- c(# ((abs(F_CbpCrCo$F_crco_logFC) <= 1) | (F_CbpCrCo$F_crco_FDR >= 0.05)) &
            ((abs(F_CbpCrCo$F_crcbp_logFC) <= 1) | (F_CbpCrCo$F_crcbp_FDR >= 0.05)) &
            ((abs(F_CbpCrCo$F_cocbp_logFC) <= 1) | (F_CbpCrCo$F_cocbp_FDR >= 0.05)))
FS_NoDifference <- F_CbpCrCo[FS_ND,]
message("In flower, number of genes with expression pattern of No Difference is: ", dim(FS_NoDifference)[1], 
        ", (", round(dim(FS_NoDifference)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
FS_Ad_1 <- c(# (F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
              ((F_CbpCrCo$F_crcbp_logFC > 1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
              ((F_CbpCrCo$F_cocbp_logFC < -1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_Ad_2 <- c(# (F_CbpCrCo$F_crco_logFC < -1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
              ((F_CbpCrCo$F_crcbp_logFC < -1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
              ((F_CbpCrCo$F_cocbp_logFC > 1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_Additivity <- F_CbpCrCo[FS_Ad_1 | FS_Ad_2,]
message("In flower, number of genes with expression pattern of Additivity is: ", dim(FS_Additivity)[1],
        ", (", round(dim(FS_Additivity)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
FS_Dom_1 <- c( #(F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
               ((abs(F_CbpCrCo$F_crcbp_logFC) <= 1) | (F_CbpCrCo$F_crcbp_FDR >= 0.05)) &
               ((F_CbpCrCo$F_cocbp_logFC < -1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_Dom_2 <- c( #F_CbpCrCo$F_crco_logFC < -1 & F_CbpCrCo$F_crco_FDR < 0.05 &
               ((abs(F_CbpCrCo$F_crcbp_logFC) <= 1) | (F_CbpCrCo$F_crcbp_FDR >= 0.05)) &
               ((F_CbpCrCo$F_cocbp_logFC > 1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_Dom_3 <- c( #(F_CbpCrCo$F_crco_logFC > 1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
               ((F_CbpCrCo$F_crcbp_logFC > 1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
               ((abs(F_CbpCrCo$F_cocbp_logFC) <= 1) | (F_CbpCrCo$F_cocbp_FDR  >= 0.05)))
FS_Dom_4 <- c( #(F_CbpCrCo$F_crco_logFC < -1) & (F_CbpCrCo$F_crco_FDR < 0.05) &
               ((F_CbpCrCo$F_crcbp_logFC < -1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
               ((abs(F_CbpCrCo$F_cocbp_logFC) <= 1) | (F_CbpCrCo$F_cocbp_FDR >= 0.05)))

FS_Dominance <- F_CbpCrCo[FS_Dom_1 | FS_Dom_2 | FS_Dom_3 | FS_Dom_4,]
message("In flower, number of genes with expression pattern of Dominance is: ", dim(FS_Dominance)[1],
        ", (", round(dim(FS_Dominance)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
FS_PTS <- c(((F_CbpCrCo$F_crcbp_logFC < -1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
             ((F_CbpCrCo$F_cocbp_logFC < -1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_PositiveTransgressive <- F_CbpCrCo[FS_PTS,]
message("In flower, number of genes with expression pattern of Positive Transgressive is: ", dim(FS_PositiveTransgressive)[1],
        ", (", round(dim(FS_PositiveTransgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
FS_NTS <- c(((F_CbpCrCo$F_crcbp_logFC > 1) & (F_CbpCrCo$F_crcbp_FDR < 0.05)) &
             ((F_CbpCrCo$F_cocbp_logFC > 1) & (F_CbpCrCo$F_cocbp_FDR < 0.05)))
FS_NegativeTransgressive <- F_CbpCrCo[FS_NTS,]
message("In flower, number of genes with expression pattern of Negative Transgressive is: ", dim(FS_NegativeTransgressive)[1],
        ", (", round(dim(FS_NegativeTransgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")
FS_Transgressive <- F_CbpCrCo[FS_PTS | FS_NTS, ]
message("In flower, number of genes with expression pattern of Negative Transgressive is: ", dim(FS_Transgressive)[1],
        ", (", round(dim(FS_Transgressive)[1]/dim(F_CbpCrCo)[1], 3) * 100, "%).")


## Leaf, Strict
# Expression pattern of No Difference: Cr = Cbp = Co
LS_ND <- c(# ((abs(L_CbpCrCo$L_crco_logFC) <= 1) | (L_CbpCrCo$L_crco_FDR >= 0.05)) &
  ((abs(L_CbpCrCo$L_crcbp_logFC) <= 1) | (L_CbpCrCo$L_crcbp_FDR >= 0.05)) &
    ((abs(L_CbpCrCo$L_cocbp_logFC) <= 1) | (L_CbpCrCo$L_cocbp_FDR >= 0.05)))
LS_NoDifference <- LS_CbpCrCo[L_ND,]
message("In leaf, number of genes with expression pattern of No Difference is: ", dim(LS_NoDifference)[1], 
        ", (", round(dim(LS_NoDifference)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
LS_Ad_1 <- c(# (L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC > 1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC < -1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_Ad_2 <- c(# (L_CbpCrCo$L_crco_logFC < -1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC < -1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC > 1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_Additivity <- L_CbpCrCo[LS_Ad_1 | LS_Ad_2,]
message("In leaf, number of genes with expression pattern of Additivity is: ", dim(LS_Additivity)[1],
        ", (", round(dim(LS_Additivity)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
LS_Dom_1 <- c( #(L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((abs(L_CbpCrCo$L_crcbp_logFC) <= 1) | (L_CbpCrCo$L_crcbp_FDR >= 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC < -1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_Dom_2 <- c( #L_CbpCrCo$L_crco_logFC < -1 & L_CbpCrCo$L_crco_FDR < 0.05 &
  ((abs(L_CbpCrCo$L_crcbp_logFC) <= 1) | (L_CbpCrCo$L_crcbp_FDR >= 0.05)) &
    ((L_CbpCrCo$L_cocbp_logFC > 1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_Dom_3 <- c( #(L_CbpCrCo$L_crco_logFC > 1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC > 1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((abs(L_CbpCrCo$L_cocbp_logFC) <= 1) | (L_CbpCrCo$L_cocbp_FDR  >= 0.05)))
LS_Dom_4 <- c( #(L_CbpCrCo$L_crco_logFC < -1) & (L_CbpCrCo$L_crco_FDR < 0.05) &
  ((L_CbpCrCo$L_crcbp_logFC < -1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
    ((abs(L_CbpCrCo$L_cocbp_logFC) <= 1) | (L_CbpCrCo$L_cocbp_FDR >= 0.05)))

LS_Dominance <- L_CbpCrCo[LS_Dom_1 | LS_Dom_2 | LS_Dom_3 | LS_Dom_4,]
message("In leaf, number of genes with expression pattern of Dominance is: ", dim(LS_Dominance)[1],
        ", (", round(dim(LS_Dominance)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
LS_PTS <- c(((L_CbpCrCo$L_crcbp_logFC < -1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
             ((L_CbpCrCo$L_cocbp_logFC < -1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_PositiveTransgressive <- L_CbpCrCo[LS_PTS,]
message("In leaf, number of genes with expression pattern of Positive Transgressive is: ", dim(LS_PositiveTransgressive)[1],
        ", (", round(dim(LS_PositiveTransgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
LS_NTS <- c(((L_CbpCrCo$L_crcbp_logFC > 1) & (L_CbpCrCo$L_crcbp_FDR < 0.05)) &
             ((L_CbpCrCo$L_cocbp_logFC > 1) & (L_CbpCrCo$L_cocbp_FDR < 0.05)))
LS_NegativeTransgressive <- L_CbpCrCo[LS_NTS,]
message("In leaf, number of genes with expression pattern of Negative Transgressive is: ", dim(LS_NegativeTransgressive)[1],
        ", (", round(dim(LS_NegativeTransgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")
LS_Transgressive <- L_CbpCrCo[LS_PTS | LS_NTS, ]
message("In leaf, number of genes with expression pattern of Negative Transgressive is: ", dim(LS_Transgressive)[1],
        ", (", round(dim(LS_Transgressive)[1]/dim(L_CbpCrCo)[1], 3) * 100, "%).")


## Root, Strict
# Expression pattern of No Difference: Cr = Cbp = Co
RS_ND <- c(# ((abs(R_CbpCrCo$R_crco_logFC) <= 1) | (R_CbpCrCo$R_crco_FDR >= 0.05)) &
  ((abs(R_CbpCrCo$R_crcbp_logFC) <= 1) | (R_CbpCrCo$R_crcbp_FDR >= 0.05)) &
    ((abs(R_CbpCrCo$R_cocbp_logFC) <= 1) | (R_CbpCrCo$R_cocbp_FDR >= 0.05)))
RS_NoDifference <- R_CbpCrCo[RS_ND,]
message("In root, number of genes with expression pattern of No Difference is: ", dim(RS_NoDifference)[1], 
        ", (", round(dim(RS_NoDifference)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Additivity: 1.) Cr > Cbp > Co, or 2.) Cr < Cbp < Co
RS_Ad_1 <- c(# (R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC > 1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC < -1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_Ad_2 <- c(# (R_CbpCrCo$R_crco_logFC < -1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC < -1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC > 1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_Additivity <- R_CbpCrCo[RS_Ad_1 | RS_Ad_2,]
message("In root, number of genes with expression pattern of Additivity is: ", dim(RS_Additivity)[1],
        ", (", round(dim(RS_Additivity)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Dominance:
# 1.) Cr = Cbp > Co; 2.) Cr = Cbp < Co; 3.) Cr > Cbp = Co; 4.) Cr < Cbp = Co
RS_Dom_1 <- c( #(R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((abs(R_CbpCrCo$R_crcbp_logFC) <= 1) | (R_CbpCrCo$R_crcbp_FDR >= 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC < -1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_Dom_2 <- c( #R_CbpCrCo$R_crco_logFC < -1 & R_CbpCrCo$R_crco_FDR < 0.05 &
  ((abs(R_CbpCrCo$R_crcbp_logFC) <= 1) | (R_CbpCrCo$R_crcbp_FDR >= 0.05)) &
    ((R_CbpCrCo$R_cocbp_logFC > 1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_Dom_3 <- c( #(R_CbpCrCo$R_crco_logFC > 1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC > 1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((abs(R_CbpCrCo$R_cocbp_logFC) <= 1) | (R_CbpCrCo$R_cocbp_FDR  >= 0.05)))
RS_Dom_4 <- c( #(R_CbpCrCo$R_crco_logFC < -1) & (R_CbpCrCo$R_crco_FDR < 0.05) &
  ((R_CbpCrCo$R_crcbp_logFC < -1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
    ((abs(R_CbpCrCo$R_cocbp_logFC) <= 1) | (R_CbpCrCo$R_cocbp_FDR >= 0.05)))

RS_Dominance <- R_CbpCrCo[RS_Dom_1 | RS_Dom_2 | RS_Dom_3 | RS_Dom_4,]
message("In root, number of genes with expression pattern of Dominance is: ", dim(RS_Dominance)[1],
        ", (", round(dim(RS_Dominance)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")

# Expression pattern of Transgressive:
# 1.) Positive Transgressive:  Cbp > Cr & Cbp > Co
RS_PTS <- c(((R_CbpCrCo$R_crcbp_logFC < -1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
             ((R_CbpCrCo$R_cocbp_logFC < -1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_PositiveTransgressive <- R_CbpCrCo[RS_PTS,]
message("In root, number of genes with expression pattern of Positive Transgressive is: ", dim(RS_PositiveTransgressive)[1],
        ", (", round(dim(RS_PositiveTransgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
# 2.) Negative Transgressive: Cbp < Cr & Cbp < Co
RS_NTS <- c(((R_CbpCrCo$R_crcbp_logFC > 1) & (R_CbpCrCo$R_crcbp_FDR < 0.05)) &
             ((R_CbpCrCo$R_cocbp_logFC > 1) & (R_CbpCrCo$R_cocbp_FDR < 0.05)))
RS_NegativeTransgressive <- R_CbpCrCo[RS_NTS,]
message("In root, number of genes with expression pattern of Negative Transgressive is: ", dim(RS_NegativeTransgressive)[1],
        ", (", round(dim(RS_NegativeTransgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")
RS_Transgressive <- R_CbpCrCo[RS_PTS | RS_NTS, ]
message("In root, number of genes with expression pattern of Negative Transgressive is: ", dim(RS_Transgressive)[1],
        ", (", round(dim(RS_Transgressive)[1]/dim(R_CbpCrCo)[1], 3) * 100, "%).")









