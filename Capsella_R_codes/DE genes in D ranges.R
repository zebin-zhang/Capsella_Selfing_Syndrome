
########    DE genes of Outcrosser and Selfer | two D indices   ########

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
Accession <- read.csv("InputData/28genomes_annot_R.csv", sep = ";")
head(Accession)

### Define populations
######################################################
CR_F <- c(Accession$populations=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$populations=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$populations=="CO" & Accession$tissue=="flower")
OC_F <- CG_F
SE_F <- c(CR_F | CO_F)

CR_L <- c(Accession$populations=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$populations=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$populations=="CO" & Accession$tissue=="leaf")
OC_L <- CG_L
SE_L <- c(CR_L | CO_L)

CR_R <- c(Accession$populations=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$populations=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$populations=="CO" & Accession$tissue=="root")
OC_R <- CG_R
SE_R <- c(CR_R | CO_R)

######################################################
raw.RC <- raw.RC[,-(37:132)] # Remove data from Cbp
head(raw.RC)
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
head(RC)

## Filter HSE
# Filter by NA number
# For each gene in each tissue in parental species, at least 5 samples contain none-zero value.
Fkeep <- rowSums(raw.RC.filter1[,c(1:12)] >= 1, na.rm = T) >= 5
Lkeep <- rowSums(raw.RC.filter1[,c(13:24)] >= 1, na.rm = T) >= 5
Rkeep <- rowSums(raw.RC.filter1[,c(25:36)] >= 1, na.rm = T) >= 5

# merge three logical vectors into one column by cbind(## & ##)
# cbind(## , ##) will combine several columns into one dataframe. 
FLRkeep <- cbind(Fkeep & Lkeep & Rkeep) 

RC <- raw.RC.filter1[FLRkeep,] 
dim(RC)
message("Number of genes BEFORE filter2: ", dim(raw.RC.filter1)[1], 
        "\nNumber of genes AFTER filter2: ",  dim(RC)[1])

tissue<-c("F","L","R") # set tissue

# Loop start
for (i in tissue){
  print(i)
  # set a sub loop for viable name
  if (i == "F"){
    subData <- "Flower"
  }else{
    if (i == "L"){
      subData <- "Leaf"
    }else{
      subData <- "Root"
    }
  }
  
  reads <- RC[, which(sub(".*_","",colnames(RC)) == i)]
  reads <- na.omit(reads)
  reads <- reads[,c(5:8, 1:4, 9:12)] # Order columns by CG CR CO
  head(reads)
  myGroup <- c(rep("Outcrosser",4), rep("Selfer", 8))
  myGroup <- factor(myGroup, levels = c("Selfer", "Outcrosser"))
  message("my group in ", i, " is: ", myGroup)
  
  myDE <- DGEList(counts = reads, group = myGroup)
  # Normalizing the data
  myDE <- calcNormFactors(myDE, method = "TMM")
  # Estimating the Dispersion
  # This is the first major step in the analysis of DGE data
  myDE <- estimateCommonDisp(myDE)
  # Extract TMM values.
  TMM.Keep <- myDE$pseudo.counts
  # Output data
  #write.table(TMM.Keep, paste("InputData/TMM_Diploids_", i, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
  
  # Create the contrast matrix
  myDE.mat <- model.matrix(~ 0 + myDE$samples$group)
  colnames(myDE.mat) <- levels(myDE$samples$group)
  message("levels of my group in ", i, " is: ", levels(myDE$samples$group)) # check the group order.
  
  # Estimate dispersion parameter for GLM
  myDE2 <- estimateGLMCommonDisp(myDE, myDE.mat)
  myDE2 <- estimateGLMTrendedDisp(myDE2, myDE.mat)
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  # I chose the method="power" after looking at the methods that are offered: bin.spline (default if number of 
  # tags is > 200), power (default otherwise), bin.loess, and spline. We have 17,396 tags, so the default is 
  # bin.spline. When I used the bin.spline method, it was no better than estimating a common dispersion so I 
  # instead used power.
  myDE2 <- estimateGLMTagwiseDisp(myDE2, myDE.mat)
  
  ## DEG analysis 
  # Model fitting
  
  fit <- glmFit(myDE2, myDE.mat)
  
    # compare (group 1 <Selfer> - group 2 <Outcrosser>) to 0: 
    # this is equivalent to comparing group 1 to group 2
  F.contrast <- makeContrasts(Selfer - Outcrosser, levels = colnames(myDE.mat))

        lrt <- glmLRT(fit, contrast = F.contrast) # group 1 to group 2, CG to CR
      vs <- "SE/OC"

    
    # Access results tables
    Table <- lrt$table
    
    # Toptags Result with FDR
    TopTags <- topTags(lrt, n = nrow(lrt$table))$table
    
    # Annotate data according to up-regulated and down-regulated
    #### ggscatter ####
    library(ggpubr)
    #install.packages("ggExtra")
    library(ggExtra)
    
    TopTags$logFDR <- -log10(TopTags$FDR)
    
    head(TopTags)
    # Loop for  TopTags$Tissue
    if (i == "F"){
      TopTags$Tissue <- "Flower"
    }else{
      if (i == "L"){
        TopTags$Tissue <- "Leaf"
      }else{
        TopTags$Tissue <- "Root"
      }
    }
    
    TopTags$DE <- "Not DE"
    TopTags[(TopTags$logFC > 0) & (TopTags$FDR <= 0.05),]$DE <- "Up" # ！not set for FC
    TopTags[(TopTags$logFC < 0) & (TopTags$FDR <= 0.05),]$DE <- "Down"
    TopTags$DE <- factor(TopTags$DE, levels = c("Up", "Down", "Not DE"))
    TopTags$Comparison <- vs
    # write data out
    # All DE result
    write.table(TopTags, paste("OutputData/DE_", i, "_SEvsOC", ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    # Extract the DE genes with FDR < 0.05
    # FC > 2
    DE <- TopTags[(abs(TopTags$logFC) > 1),]
    write.table(DE, paste("OutputData/DE_", i, "_SEvsOC", "_FC2.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    # FDR ≤ 0.05
    DE <- TopTags[(TopTags$FDR <= 0.05),]
    write.table(DE, paste("OutputData/DE_", i, "_SEvsOC", "_FDR.05.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    assign(paste0("DE_",i,"_SEvsOC"), DE)
    # FC > 2 & FDR ≤ 0.05
    DE <- TopTags[(abs(TopTags$logFC) > 1) & (TopTags$FDR <= 0.05),]
    write.table(DE, paste("OutputData/DE_", i, "_SEvsOC", "_FC2FDR.05.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    #assign(paste0("DE_",i,"_",com), DE)
    # FC > 2 & FDR <= 0.01
    DE <- TopTags[(abs(TopTags$logFC) > 1) & (TopTags$FDR <= 0.01),]
    write.table(DE, paste("OutputData/DE_", i, "_SEvsOC", "_FC2FDR.01.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    
    ### Message out 
    count_R <- dim(TopTags[TopTags$DE == "Up",])[1]
    count_L <- dim(TopTags[TopTags$DE == "Down",])[1]
    median_R <- round(median(TopTags[TopTags$DE == "Up",]$logFC), 2)
    median_L <- round(median(TopTags[TopTags$DE == "Down",]$logFC), 2)
    percentage_R <- paste0(round((dim(TopTags[TopTags$DE == "Up",])[1] / dim(TopTags)[1])*100, 2), "%")
    percentage_L <- paste0(round((dim(TopTags[TopTags$DE == "Down",])[1] / dim(TopTags)[1])*100, 2), "%")
    
    message("Numbers of Up regulated genes of SE vs OC in ", i, " is: ", dim(TopTags[TopTags$DE == "Up",])[1])
    message("Numbers of Down regulated genes of SE vs OC in ", i, " is: ", dim(TopTags[TopTags$DE == "Down",])[1])
    message("Percentage of Up regulated genes of SE vs OC in ", i, " is: ", (dim(TopTags[TopTags$DE == "Up",])[1] / dim(TopTags)[1])*100)
    message("Percentage of Down regulated genes of SE vs OC in ", i, " is: ", (dim(TopTags[TopTags$DE == "Down",])[1] / dim(TopTags)[1])*100)
  
}

###############################################################
#                  D indices: Dco & Dcr                       #
###############################################################

library(stringr)
library(matrixStats)
library(ggplot2)
library(cowplot)

for (tissue in c("F", "L", "R")) {
  Reads <- read.table(paste("InputData/TMM_Diploids_", tissue, ".txt", sep = ""), header = T, sep = "\t")
  head(Reads)
  CGColumns <- str_detect(names(Reads), "CG") # detect which columns contain the "CG" 
  CRColumns <- str_detect(names(Reads), "CR") # detect which columns contain the "CG" 
  COColumns <- str_detect(names(Reads), "CO") # detect which columns contain the "CG" 
  
  Reads$MeanCG <- rowMeans(Reads[, CGColumns], na.rm = T)
  Reads$MeanCR <- rowMeans(Reads[, CRColumns], na.rm = T)
  Reads$MeanCO <- rowMeans(Reads[, COColumns], na.rm = T)
  
  Reads$D.RO <- abs(Reads$MeanCR - Reads$MeanCO) # Distance from CO to CR
  Reads$D.RG <- abs(Reads$MeanCR - Reads$MeanCG) # Distance from CO to CG
  Reads$D.OR <- abs(Reads$MeanCO - Reads$MeanCR) # Distance from CO to CR
  Reads$D.OG <- abs(Reads$MeanCO - Reads$MeanCG) # Distance from CO to CG
  Reads$D.R_Max <- rowMaxs(as.matrix(Reads[,c("D.RO","D.RG")]), na.rm = TRUE)
  Reads$D.O_Max <- rowMaxs(as.matrix(Reads[,c("D.OR","D.OG")]), na.rm = TRUE)
  
  Reads$D.R <- (Reads$D.RO - Reads$D.RG) / Reads$D.R_Max
  Reads$D.O <- (Reads$D.OR - Reads$D.OG) / Reads$D.O_Max
  
  assign(paste0("Reads_", tissue), Reads)
}

# combine DE genes with D values
head(DE_F_SEvsOC)
head(Reads_F)
D_F <-  Reads_F[,c("D.O","D.R")]
DE_F_SEvsOC_D <- merge(DE_F_SEvsOC, Reads_F, by = 0) # merge by row names
DE_L_SEvsOC_D <- merge(DE_L_SEvsOC, Reads_L, by = 0) # merge by row names
DE_R_SEvsOC_D <- merge(DE_R_SEvsOC, Reads_R, by = 0) # merge by row names
head(DE_F_SEvsOC_D)
dim(DE_F_SEvsOC)
dim(DE_F_SEvsOC_D) 
str(DE_F_SEvsOC_D)
library(ggplot2)
library(scales)
dim(DE_F_SEvsOC_D)
p1 <- ggplot(DE_F_SEvsOC_D, aes(D.O, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#bb99b7","#ecc8c9")) +
  xlim(-1,1) +
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=1577) in flower", x = expression(D["CO"]), y = "Percentage") + 
  theme_classic()

p2 <- ggplot(DE_F_SEvsOC_D, aes(D.R, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#e48826","#c6a78f")) +
  xlim(-1,1)+
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=1577) in flower",x = expression(D["CR"]), y = "Percentage") + 
  theme_classic()

dim(DE_L_SEvsOC_D)
p3 <- ggplot(DE_L_SEvsOC_D, aes(D.O, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#bb99b7","#ecc8c9")) +
  xlim(-1,1)+
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=213) in leaf",x = expression(D["CO"]), y = "Percentage") + 
  theme_classic()

p4 <- ggplot(DE_L_SEvsOC_D, aes(D.R, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#e48826","#c6a78f")) +
  xlim(-1,1)+
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=213) in leaf",x = expression(D["CR"]), y = "Percentage") + 
  theme_classic()

dim(DE_R_SEvsOC_D)
p5 <- ggplot(DE_R_SEvsOC_D, aes(D.O, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#bb99b7","#ecc8c9")) +
  xlim(-1,1)+
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=377) in root",x = expression(D["CO"]), y = "Percentage") + 
  theme_classic()

p6 <- ggplot(DE_R_SEvsOC_D, aes(D.R, fill = DE)) +
  geom_histogram(bins = 20, position = "fill", color = "black") +
  scale_fill_manual(values = c("#e48826","#c6a78f")) +
  xlim(-1,1)+
  scale_y_continuous(labels = percent) +
  labs(title = "Selfers/Outcrossers DE genes (n=377) in root",x = expression(D["CR"]), y = "Percentage") + 
  theme_classic()

library(ggpubr)
library("grid")
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,2)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p3, vp = vplayout(2,1))
print(p4, vp = vplayout(2,2))
print(p5, vp = vplayout(3,1))
print(p6, vp = vplayout(3,2))




