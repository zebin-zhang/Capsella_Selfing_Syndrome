###############################################################
#          Similarity and convergence indices                 #
###############################################################

# To quantify the similarity between each subgenome expression level and the expression level in the parental 
# species, we developed a similarity index (S). For each transcript i and each subgenome j in {CbpCg, CbpCo}, 
# S was computed as the subgenome relative expression deviation from the mean expression level in the parental
# species, u(i) = (E(ico) + E(icg)) / 2:
# S(ij) = (E(ij) - u(i)) / u(i)

## Set work directory and load files
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

TMM <- read.table("InputData/Capsella_ReadsCount_TMM_Filtered.txt", header = T, sep = "\t")
head(TMM)
dim(TMM)
TMM[TMM == 0] <- NA

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
######################################################

# Similarity index

# Create a new dataframe for Flower
S_Flower <- data.frame(Gene=rownames(TMM))
S_Flower$mean_CG_F <- rowMeans(TMM[,CG_F], na.rm = T) # Average value of CG expression in Flower
S_Flower$mean_CR_F <- rowMeans(TMM[,CR_F], na.rm = T)
S_Flower$mean_CO_F <- rowMeans(TMM[,CO_F], na.rm = T)

S_Flower$mean_Cbp_Cg_F <- rowMeans(TMM[,Cbp_Cg_F], na.rm = T)
S_Flower$mean_Cbp_Co_F <- rowMeans(TMM[,Cbp_Co_F], na.rm = T)

S_Flower$mean_ASI_Cg_F <- rowMeans(TMM[,ASI_Cg_F], na.rm = T)
S_Flower$mean_EUR_Cg_F <- rowMeans(TMM[,EUR_Cg_F], na.rm = T)
S_Flower$mean_ME_Cg_F <- rowMeans(TMM[,ME_Cg_F], na.rm = T)
S_Flower$mean_CASI_Cg_F <- rowMeans(TMM[,CASI_Cg_F], na.rm = T)

S_Flower$mean_ASI_Co_F <- rowMeans(TMM[,ASI_Co_F], na.rm = T)
S_Flower$mean_EUR_Co_F <- rowMeans(TMM[,EUR_Co_F], na.rm = T)
S_Flower$mean_ME_Co_F <- rowMeans(TMM[,ME_Co_F], na.rm = T)
S_Flower$mean_CASI_Co_F <- rowMeans(TMM[,CASI_Co_F], na.rm = T)

head(TMM[,CG_F])
head(TMM[,EUR_Cg_F])
head(S_Flower)
# write.table(S_Flower, file = "OutputData/meanTMM_of_phased_Cbp_in_F.txt", quote = F, sep = "\t", row.names = T, col.names = T)
dim(S_Flower)

# Average expression level in parents of CR & CO.
S_Flower$Ave_CRCO <- rowMeans(S_Flower[,c("mean_CR_F","mean_CO_F")])

# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
S_Flower$Scrco_Cbp_Cg <- (S_Flower$mean_Cbp_Cg_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO # Subgenome Cg 

S_Flower$Scrco_ASI_Cg <- (S_Flower$mean_ASI_Cg_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_EUR_Cg <- (S_Flower$mean_EUR_Cg_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_ME_Cg <- (S_Flower$mean_ME_Cg_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_CASI_Cg <- (S_Flower$mean_CASI_Cg_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO

S_Flower$Scrco_Cbp_Co <- (S_Flower$mean_Cbp_Co_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO # Subgenome Co

S_Flower$Scrco_ASI_Co <- (S_Flower$mean_ASI_Co_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_EUR_Co <- (S_Flower$mean_EUR_Co_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_ME_Co <- (S_Flower$mean_ME_Co_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO
S_Flower$Scrco_CASI_Co <- (S_Flower$mean_CASI_Co_F - S_Flower$Ave_CRCO) / S_Flower$Ave_CRCO

head(S_Flower)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_Cbp_Cg <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_Cbp_Cg * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ASI_Cg <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ASI_Cg * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_EUR_Cg <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_EUR_Cg * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ME_Cg <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ME_Cg * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_CASI_Cg <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_CASI_Cg * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_Cbp_Co <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_Cbp_Co * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ASI_Co <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ASI_Co * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_EUR_Co <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_EUR_Co * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ME_Co <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_ME_Co * (-1)
S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_CASI_Co <- S_Flower[(S_Flower$mean_CR_F > S_Flower$mean_CO_F),]$Scrco_CASI_Co * (-1)

# Calculate delta S by comparing the absolute values of S_Cbp_Cr and S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
S_Flower$delta_Scrco_Cbp <- abs(S_Flower$Scrco_Cbp_Co) - abs(S_Flower$Scrco_Cbp_Cg)
S_Flower$delta_Scrco_ASI <- abs(S_Flower$Scrco_ASI_Co) - abs(S_Flower$Scrco_ASI_Cg)
S_Flower$delta_Scrco_EUR <- abs(S_Flower$Scrco_EUR_Co) - abs(S_Flower$Scrco_EUR_Cg)
S_Flower$delta_Scrco_ME <- abs(S_Flower$Scrco_ME_Co) - abs(S_Flower$Scrco_ME_Cg)
S_Flower$delta_Scrco_CASI <- abs(S_Flower$Scrco_CASI_Co) - abs(S_Flower$Scrco_CASI_Cg)
head(S_Flower)

# load genes significantly expressed between Cr Co in flower.
DE.F.crco <- read.table("OutputData/DE_CrCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.crco)
dim(DE.F.crco)

# keep genes with FDR ≤ 0.01
DE.F.crco.keep <- S_Flower$Gene %in% row.names(DE.F.crco)
S_Flower$Gene %in% row.names(DE.F.crco)

Scrco_Flower_FDR <- S_Flower[DE.F.crco.keep,]
head(Scrco_Flower_FDR)
dim(Scrco_Flower_FDR)
##################################

# Create a new dataframe for Leaf
S_Leaf <- data.frame(Gene=rownames(TMM))
S_Leaf$mean_CG_L <- rowMeans(TMM[,CG_L], na.rm = T) # Average value of CG expression in Leaf
S_Leaf$mean_CR_L <- rowMeans(TMM[,CR_L], na.rm = T)
S_Leaf$mean_CO_L <- rowMeans(TMM[,CO_L], na.rm = T)

S_Leaf$mean_Cbp_Cg_L <- rowMeans(TMM[,Cbp_Cg_L], na.rm = T)
S_Leaf$mean_Cbp_Co_L <- rowMeans(TMM[,Cbp_Co_L], na.rm = T)

S_Leaf$mean_ASI_Cg_L <- rowMeans(TMM[,ASI_Cg_L], na.rm = T)
S_Leaf$mean_EUR_Cg_L <- rowMeans(TMM[,EUR_Cg_L], na.rm = T)
S_Leaf$mean_ME_Cg_L <- rowMeans(TMM[,ME_Cg_L], na.rm = T)
S_Leaf$mean_CASI_Cg_L <- rowMeans(TMM[,CASI_Cg_L], na.rm = T)

S_Leaf$mean_ASI_Co_L <- rowMeans(TMM[,ASI_Co_L], na.rm = T)
S_Leaf$mean_EUR_Co_L <- rowMeans(TMM[,EUR_Co_L], na.rm = T)
S_Leaf$mean_ME_Co_L <- rowMeans(TMM[,ME_Co_L], na.rm = T)
S_Leaf$mean_CASI_Co_L <- rowMeans(TMM[,CASI_Co_L], na.rm = T)

head(TMM[,CG_L])
head(TMM[,EUR_Cg_L])
head(S_Leaf)
# write.table(S_Leaf, file = "OutputData/meanTMM_of_phased_Cbp_in_L.txt", quote = F, sep = "\t", row.names = T, col.names = T)

# Average expression level in parents of CR & CO.
S_Leaf$Ave_CRCO <- rowMeans(S_Leaf[,c("mean_CR_L","mean_CO_L")])
S_Leaf$Ave_CGCO <- rowMeans(S_Leaf[,c("mean_CG_L","mean_CO_L")])
S_Leaf$Ave_CGCR <- rowMeans(S_Leaf[,c("mean_CG_L","mean_CR_L")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
S_Leaf$Scrco_Cbp_Cg <- (S_Leaf$mean_Cbp_Cg_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO # Subgenome Cg 

S_Leaf$Scrco_ASI_Cg <- (S_Leaf$mean_ASI_Cg_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_EUR_Cg <- (S_Leaf$mean_EUR_Cg_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_ME_Cg <- (S_Leaf$mean_ME_Cg_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_CASI_Cg <- (S_Leaf$mean_CASI_Cg_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO

S_Leaf$Scrco_Cbp_Co <- (S_Leaf$mean_Cbp_Co_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO # Subgenome Co

S_Leaf$Scrco_ASI_Co <- (S_Leaf$mean_ASI_Co_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_EUR_Co <- (S_Leaf$mean_EUR_Co_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_ME_Co <- (S_Leaf$mean_ME_Co_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO
S_Leaf$Scrco_CASI_Co <- (S_Leaf$mean_CASI_Co_L - S_Leaf$Ave_CRCO) / S_Leaf$Ave_CRCO

head(S_Leaf)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.
for (i in c(1:dim(S_Leaf)[1])) {
  if (S_Leaf[i,]$mean_CR_L > S_Leaf[i,]$mean_CO_L) {
    S_Leaf[i,]$Scrco_Cbp_Cg <- S_Leaf[i,]$Scrco_Cbp_Cg * (-1)
    S_Leaf[i,]$Scrco_ASI_Cg <- S_Leaf[i,]$Scrco_ASI_Cg * (-1)
    S_Leaf[i,]$Scrco_EUR_Cg <- S_Leaf[i,]$Scrco_EUR_Cg * (-1)
    S_Leaf[i,]$Scrco_ME_Cg <- S_Leaf[i,]$Scrco_ME_Cg * (-1)
    S_Leaf[i,]$Scrco_CASI_Cg <- S_Leaf[i,]$Scrco_CASI_Cg * (-1)
    S_Leaf[i,]$Scrco_Cbp_Co <- S_Leaf[i,]$Scrco_Cbp_Co * (-1)
    S_Leaf[i,]$Scrco_ASI_Co <- S_Leaf[i,]$Scrco_ASI_Co * (-1)
    S_Leaf[i,]$Scrco_EUR_Co <- S_Leaf[i,]$Scrco_EUR_Co * (-1)
    S_Leaf[i,]$Scrco_ME_Co <- S_Leaf[i,]$Scrco_ME_Co * (-1)
    S_Leaf[i,]$Scrco_CASI_Co <- S_Leaf[i,]$Scrco_CASI_Co * (-1)
  }
}

# Calculate delta S by comparing the absolute values of S_Cbp_Cr and S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
S_Leaf$delta_Scrco_Cbp <- abs(S_Leaf$Scrco_Cbp_Co) - abs(S_Leaf$Scrco_Cbp_Cg)
S_Leaf$delta_Scrco_ASI <- abs(S_Leaf$Scrco_ASI_Co) - abs(S_Leaf$Scrco_ASI_Cg)
S_Leaf$delta_Scrco_EUR <- abs(S_Leaf$Scrco_EUR_Co) - abs(S_Leaf$Scrco_EUR_Cg)
S_Leaf$delta_Scrco_ME <- abs(S_Leaf$Scrco_ME_Co) - abs(S_Leaf$Scrco_ME_Cg)
S_Leaf$delta_Scrco_CASI <- abs(S_Leaf$Scrco_CASI_Co) - abs(S_Leaf$Scrco_CASI_Cg)
head(S_Leaf)

# load genes significantly expressed between Cr Co in flower.
DE.L.crco <- read.table("OutputData/DE_CrCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.crco)
dim(DE.L.crco)

# keep genes with FDR ≤ 0.01
DE.L.crco.keep <- S_Leaf$Gene %in% row.names(DE.L.crco)
S_Leaf$Gene %in% row.names(DE.L.crco)

Scrco_Leaf_FDR <- S_Leaf[DE.L.crco.keep,]
head(Scrco_Leaf_FDR)
dim(Scrco_Leaf_FDR)

#############################

# Create a new dataframe for Root
S_Root <- data.frame(Gene=rownames(TMM))
S_Root$mean_CG_R <- rowMeans(TMM[,CG_R], na.rm = T) # Average value of CG expression in Root
S_Root$mean_CR_R <- rowMeans(TMM[,CR_R], na.rm = T)
S_Root$mean_CO_R <- rowMeans(TMM[,CO_R], na.rm = T)

S_Root$mean_Cbp_Cg_R <- rowMeans(TMM[,Cbp_Cg_R], na.rm = T)
S_Root$mean_Cbp_Co_R <- rowMeans(TMM[,Cbp_Co_R], na.rm = T)

S_Root$mean_ASI_Cg_R <- rowMeans(TMM[,ASI_Cg_R], na.rm = T)
S_Root$mean_EUR_Cg_R <- rowMeans(TMM[,EUR_Cg_R], na.rm = T)
S_Root$mean_ME_Cg_R <- rowMeans(TMM[,ME_Cg_R], na.rm = T)
S_Root$mean_CASI_Cg_R <- rowMeans(TMM[,CASI_Cg_R], na.rm = T)

S_Root$mean_ASI_Co_R <- rowMeans(TMM[,ASI_Co_R], na.rm = T)
S_Root$mean_EUR_Co_R <- rowMeans(TMM[,EUR_Co_R], na.rm = T)
S_Root$mean_ME_Co_R <- rowMeans(TMM[,ME_Co_R], na.rm = T)
S_Root$mean_CASI_Co_R <- rowMeans(TMM[,CASI_Co_R], na.rm = T)

head(TMM[,CG_R])
head(TMM[,EUR_Cg_R])
head(S_Root)
# write.table(S_Root, file = "OutputData/meanTMM_of_phased_Cbp_in_R.txt", quote = F, sep = "\t", row.names = T, col.names = T)

# Average expression level in parents of CR & CO.
S_Root$Ave_CRCO <- rowMeans(S_Root[,c("mean_CR_R","mean_CO_R")])
S_Root$Ave_CGCO <- rowMeans(S_Root[,c("mean_CG_R","mean_CO_R")])
S_Root$Ave_CGCR <- rowMeans(S_Root[,c("mean_CG_R","mean_CR_R")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
S_Root$Scrco_Cbp_Cg <- (S_Root$mean_Cbp_Cg_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO # Subgenome Cg 

S_Root$Scrco_ASI_Cg <- (S_Root$mean_ASI_Cg_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_EUR_Cg <- (S_Root$mean_EUR_Cg_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_ME_Cg <- (S_Root$mean_ME_Cg_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_CASI_Cg <- (S_Root$mean_CASI_Cg_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO

S_Root$Scrco_Cbp_Co <- (S_Root$mean_Cbp_Co_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO # Subgenome Co

S_Root$Scrco_ASI_Co <- (S_Root$mean_ASI_Co_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_EUR_Co <- (S_Root$mean_EUR_Co_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_ME_Co <- (S_Root$mean_ME_Co_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO
S_Root$Scrco_CASI_Co <- (S_Root$mean_CASI_Co_R - S_Root$Ave_CRCO) / S_Root$Ave_CRCO

head(S_Root)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.
orientate <- (S_Root$mean_CR_R > S_Root$mean_CO_R)

    S_Root[orientate,]$Scrco_Cbp_Cg <- S_Root[orientate,]$Scrco_Cbp_Cg * (-1)
    S_Root[orientate,]$Scrco_ASI_Cg <- S_Root[orientate,]$Scrco_ASI_Cg * (-1)
    S_Root[orientate,]$Scrco_EUR_Cg <- S_Root[orientate,]$Scrco_EUR_Cg * (-1)
    S_Root[orientate,]$Scrco_ME_Cg <- S_Root[orientate,]$Scrco_ME_Cg * (-1)
    S_Root[orientate,]$Scrco_CASI_Cg <- S_Root[orientate,]$Scrco_CASI_Cg * (-1)
    S_Root[orientate,]$Scrco_Cbp_Co <- S_Root[orientate,]$Scrco_Cbp_Co * (-1)
    S_Root[orientate,]$Scrco_ASI_Co <- S_Root[orientate,]$Scrco_ASI_Co * (-1)
    S_Root[orientate,]$Scrco_EUR_Co <- S_Root[orientate,]$Scrco_EUR_Co * (-1)
    S_Root[orientate,]$Scrco_ME_Co <- S_Root[orientate,]$Scrco_ME_Co * (-1)
    S_Root[orientate,]$Scrco_CASI_Co <- S_Root[orientate,]$Scrco_CASI_Co * (-1)

# Calculate delta S by comparing the absolute values of S_Cbp_Cr and S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
S_Root$delta_Scrco_Cbp <- abs(S_Root$Scrco_Cbp_Co) - abs(S_Root$Scrco_Cbp_Cg)
S_Root$delta_Scrco_ASI <- abs(S_Root$Scrco_ASI_Co) - abs(S_Root$Scrco_ASI_Cg)
S_Root$delta_Scrco_EUR <- abs(S_Root$Scrco_EUR_Co) - abs(S_Root$Scrco_EUR_Cg)
S_Root$delta_Scrco_ME <- abs(S_Root$Scrco_ME_Co) - abs(S_Root$Scrco_ME_Cg)
S_Root$delta_Scrco_CASI <- abs(S_Root$Scrco_CASI_Co) - abs(S_Root$Scrco_CASI_Cg)
head(S_Root)


# load genes significantly expressed between Cr Co in flower.
DE.R.crco <- read.table("OutputData/DE_CrCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.crco)
dim(DE.R.crco)


# keep genes with FDR ≤ 0.01
DE.R.crco.keep <- S_Root$Gene %in% row.names(DE.R.crco)
S_Root$Gene %in% row.names(DE.R.crco)

Scrco_Root_FDR <- S_Root[DE.R.crco.keep,]
head(Scrco_Root_FDR)
dim(Scrco_Root_FDR)

######################################################

# Convergence index
# for each transcript i, we used a convergence index Ci = ∆parents - ∆x / max(∆parents, ∆x)
# ∆x stands for either ΔCo, ΔCg or Δsub
# Ci thus ranges from -1 to 1, with positive values indicating more similar expression between the subgenomes
# of C. bursa-pastoris than between parental species, and negative values indicating increased diffeRences
# between subgenomes; the closer Ci to 0, the more similar are the expression patterns to parental species.

library(MatrixGenerics)

# Flower
head(Scrco_Flower_FDR)

# Within dataframe of differently expressed between parental species. CR CO
# Absolute differences between parents in expression.
Scrco_Flower_FDR$C_delta_Parents <- abs(Scrco_Flower_FDR$mean_CR_F - Scrco_Flower_FDR$mean_CO_F)

# Absolute differences between Subgenomes in expression.
Scrco_Flower_FDR$C_delta_Cbp <- abs(Scrco_Flower_FDR$mean_Cbp_Cg_F - Scrco_Flower_FDR$mean_Cbp_Co_F) # Cbp all
Scrco_Flower_FDR$C_delta_ASI <- abs(Scrco_Flower_FDR$mean_ASI_Cg_F - Scrco_Flower_FDR$mean_ASI_Co_F) # Cbp-ASI
Scrco_Flower_FDR$C_delta_EUR <- abs(Scrco_Flower_FDR$mean_EUR_Cg_F - Scrco_Flower_FDR$mean_EUR_Co_F) # Cbp-EUR
Scrco_Flower_FDR$C_delta_ME <- abs(Scrco_Flower_FDR$mean_ME_Cg_F - Scrco_Flower_FDR$mean_ME_Co_F) # Cbp-ME
Scrco_Flower_FDR$C_delta_CASI <- abs(Scrco_Flower_FDR$mean_CASI_Cg_F - Scrco_Flower_FDR$mean_CASI_Co_F) # Cbp-CASI

# each subgenome and the opposite parental species
Scrco_Flower_FDR$C_delta_Cbp_cg <- abs(Scrco_Flower_FDR$mean_Cbp_Cg_F - Scrco_Flower_FDR$mean_CO_F) # Cbp all
Scrco_Flower_FDR$C_delta_Cbp_co <- abs(Scrco_Flower_FDR$mean_Cbp_Co_F - Scrco_Flower_FDR$mean_CR_F)

Scrco_Flower_FDR$C_delta_ASI_cg <- abs(Scrco_Flower_FDR$mean_ASI_Cg_F - Scrco_Flower_FDR$mean_CO_F) # ASI
Scrco_Flower_FDR$C_delta_ASI_co <- abs(Scrco_Flower_FDR$mean_ASI_Co_F - Scrco_Flower_FDR$mean_CR_F)

Scrco_Flower_FDR$C_delta_EUR_cg <- abs(Scrco_Flower_FDR$mean_EUR_Cg_F - Scrco_Flower_FDR$mean_CO_F) # EUR
Scrco_Flower_FDR$C_delta_EUR_co <- abs(Scrco_Flower_FDR$mean_EUR_Co_F - Scrco_Flower_FDR$mean_CR_F)

Scrco_Flower_FDR$C_delta_ME_cg <- abs(Scrco_Flower_FDR$mean_ME_Cg_F - Scrco_Flower_FDR$mean_CO_F) # ME
Scrco_Flower_FDR$C_delta_ME_co <- abs(Scrco_Flower_FDR$mean_ME_Co_F - Scrco_Flower_FDR$mean_CR_F)

Scrco_Flower_FDR$C_delta_CASI_cg <- abs(Scrco_Flower_FDR$mean_CASI_Cg_F - Scrco_Flower_FDR$mean_CO_F) # CASI
Scrco_Flower_FDR$C_delta_CASI_co <- abs(Scrco_Flower_FDR$mean_CASI_Co_F - Scrco_Flower_FDR$mean_CR_F)

head(Scrco_Flower_FDR)

# define columns with empty value.
Scrco_Flower_FDR$C.Cbp <- NA
Scrco_Flower_FDR$C.Cbp_cg <- NA
Scrco_Flower_FDR$C.Cbp_co <- NA

Scrco_Flower_FDR$C.ASI <- NA
Scrco_Flower_FDR$C.ASI_cg <- NA
Scrco_Flower_FDR$C.ASI_co <- NA

Scrco_Flower_FDR$C.EUR <- NA
Scrco_Flower_FDR$C.EUR_cg <- NA
Scrco_Flower_FDR$C.EUR_co <- NA

Scrco_Flower_FDR$C.ME <- NA
Scrco_Flower_FDR$C.ME_cg <- NA
Scrco_Flower_FDR$C.ME_co <- NA

Scrco_Flower_FDR$C.CASI <- NA
Scrco_Flower_FDR$C.CASI_cg <- NA
Scrco_Flower_FDR$C.CASI_co <- NA

library("matrixStats")

for(i in c(1:dim(Scrco_Flower_FDR)[1])){
  Scrco_Flower_FDR[i,]$C.Cbp <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_Cbp) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_Cbp")]))
  Scrco_Flower_FDR[i,]$C.Cbp_cg <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_Cbp_cg) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_Cbp_cg")]))
  Scrco_Flower_FDR[i,]$C.Cbp_co <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_Cbp_co) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_Cbp_co")]))
  
  Scrco_Flower_FDR[i,]$C.ASI <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ASI) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ASI")]))
  Scrco_Flower_FDR[i,]$C.ASI_cg <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ASI_cg) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ASI_cg")]))
  Scrco_Flower_FDR[i,]$C.ASI_co <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ASI_co) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ASI_co")]))
  
  Scrco_Flower_FDR[i,]$C.EUR <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_EUR) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_EUR")]))
  Scrco_Flower_FDR[i,]$C.EUR_cg <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_EUR_cg) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_EUR_cg")]))
  Scrco_Flower_FDR[i,]$C.EUR_co <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_EUR_co) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_EUR_co")]))
  
  Scrco_Flower_FDR[i,]$C.ME <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ME) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ME")]))
  Scrco_Flower_FDR[i,]$C.ME_cg <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ME_cg) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ME_cg")]))
  Scrco_Flower_FDR[i,]$C.ME_co <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_ME_co) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_ME_co")]))
  
  Scrco_Flower_FDR[i,]$C.CASI <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_CASI) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_CASI")]))
  Scrco_Flower_FDR[i,]$C.CASI_cg <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_CASI_cg) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_CASI_cg")]))
  Scrco_Flower_FDR[i,]$C.CASI_co <- (Scrco_Flower_FDR[i,]$C_delta_Parents - Scrco_Flower_FDR[i,]$C_delta_CASI_co) / rowMaxs(as.matrix(Scrco_Flower_FDR[i, c("C_delta_Parents","C_delta_CASI_co")]))
}

head(Scrco_Flower_FDR)

# Leaf
head(Scrco_Leaf_FDR)

# Within dataframe of differently expressed between parental species. CR CO
# Absolute differences between parents in expression.
Scrco_Leaf_FDR$C_delta_Parents <- abs(Scrco_Leaf_FDR$mean_CR_L - Scrco_Leaf_FDR$mean_CO_L)

# Absolute differences between Subgenomes in expression.
Scrco_Leaf_FDR$C_delta_Cbp <- abs(Scrco_Leaf_FDR$mean_Cbp_Cg_L - Scrco_Leaf_FDR$mean_Cbp_Co_L) # Cbp all
Scrco_Leaf_FDR$C_delta_ASI <- abs(Scrco_Leaf_FDR$mean_ASI_Cg_L - Scrco_Leaf_FDR$mean_ASI_Co_L) # Cbp-ASI
Scrco_Leaf_FDR$C_delta_EUR <- abs(Scrco_Leaf_FDR$mean_EUR_Cg_L - Scrco_Leaf_FDR$mean_EUR_Co_L) # Cbp-EUR
Scrco_Leaf_FDR$C_delta_ME <- abs(Scrco_Leaf_FDR$mean_ME_Cg_L - Scrco_Leaf_FDR$mean_ME_Co_L) # Cbp-ME
Scrco_Leaf_FDR$C_delta_CASI <- abs(Scrco_Leaf_FDR$mean_CASI_Cg_L - Scrco_Leaf_FDR$mean_CASI_Co_L) # Cbp-CASI

# each subgenome and the opposite parental species
Scrco_Leaf_FDR$C_delta_Cbp_cg <- abs(Scrco_Leaf_FDR$mean_Cbp_Cg_L - Scrco_Leaf_FDR$mean_CO_L) # Cbp all
Scrco_Leaf_FDR$C_delta_Cbp_co <- abs(Scrco_Leaf_FDR$mean_Cbp_Co_L - Scrco_Leaf_FDR$mean_CR_L)

Scrco_Leaf_FDR$C_delta_ASI_cg <- abs(Scrco_Leaf_FDR$mean_ASI_Cg_L - Scrco_Leaf_FDR$mean_CO_L) # ASI
Scrco_Leaf_FDR$C_delta_ASI_co <- abs(Scrco_Leaf_FDR$mean_ASI_Co_L - Scrco_Leaf_FDR$mean_CR_L)

Scrco_Leaf_FDR$C_delta_EUR_cg <- abs(Scrco_Leaf_FDR$mean_EUR_Cg_L - Scrco_Leaf_FDR$mean_CO_L) # EUR
Scrco_Leaf_FDR$C_delta_EUR_co <- abs(Scrco_Leaf_FDR$mean_EUR_Co_L - Scrco_Leaf_FDR$mean_CR_L)

Scrco_Leaf_FDR$C_delta_ME_cg <- abs(Scrco_Leaf_FDR$mean_ME_Cg_L - Scrco_Leaf_FDR$mean_CO_L) # ME
Scrco_Leaf_FDR$C_delta_ME_co <- abs(Scrco_Leaf_FDR$mean_ME_Co_L - Scrco_Leaf_FDR$mean_CR_L)

Scrco_Leaf_FDR$C_delta_CASI_cg <- abs(Scrco_Leaf_FDR$mean_CASI_Cg_L - Scrco_Leaf_FDR$mean_CO_L) # CASI
Scrco_Leaf_FDR$C_delta_CASI_co <- abs(Scrco_Leaf_FDR$mean_CASI_Co_L - Scrco_Leaf_FDR$mean_CR_L)

head(Scrco_Leaf_FDR)

# define columns with empty value.
Scrco_Leaf_FDR$C.Cbp <- NA
Scrco_Leaf_FDR$C.Cbp_cg <- NA
Scrco_Leaf_FDR$C.Cbp_co <- NA

Scrco_Leaf_FDR$C.ASI <- NA
Scrco_Leaf_FDR$C.ASI_cg <- NA
Scrco_Leaf_FDR$C.ASI_co <- NA

Scrco_Leaf_FDR$C.EUR <- NA
Scrco_Leaf_FDR$C.EUR_cg<- NA
Scrco_Leaf_FDR$C.EUR_co <- NA

Scrco_Leaf_FDR$C.ME <- NA
Scrco_Leaf_FDR$C.ME_cg <- NA
Scrco_Leaf_FDR$C.ME_co <- NA

Scrco_Leaf_FDR$C.CASI <- NA
Scrco_Leaf_FDR$C.CASI_cg <- NA
Scrco_Leaf_FDR$C.CASI_co <- NA

for(i in c(1:dim(Scrco_Leaf_FDR)[1])){
  Scrco_Leaf_FDR[i,]$C.Cbp <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_Cbp) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_Cbp")]))
  Scrco_Leaf_FDR[i,]$C.Cbp_cg <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_Cbp_cg) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_Cbp_cg")]))
  Scrco_Leaf_FDR[i,]$C.Cbp_co <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_Cbp_co) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_Cbp_co")]))
  
  Scrco_Leaf_FDR[i,]$C.ASI <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ASI) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ASI")]))
  Scrco_Leaf_FDR[i,]$C.ASI_cg <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ASI_cg) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ASI_cg")]))
  Scrco_Leaf_FDR[i,]$C.ASI_co <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ASI_co) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ASI_co")]))
  
  Scrco_Leaf_FDR[i,]$C.EUR <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_EUR) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_EUR")]))
  Scrco_Leaf_FDR[i,]$C.EUR_cg <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_EUR_cg) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_EUR_cg")]))
  Scrco_Leaf_FDR[i,]$C.EUR_co <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_EUR_co) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_EUR_co")]))
  
  Scrco_Leaf_FDR[i,]$C.ME <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ME) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ME")]))
  Scrco_Leaf_FDR[i,]$C.ME_cg <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ME_cg) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ME_cg")]))
  Scrco_Leaf_FDR[i,]$C.ME_co <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_ME_co) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_ME_co")]))
  
  Scrco_Leaf_FDR[i,]$C.CASI <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_CASI) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_CASI")]))
  Scrco_Leaf_FDR[i,]$C.CASI_cg <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_CASI_cg) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_CASI_cg")]))
  Scrco_Leaf_FDR[i,]$C.CASI_co <- (Scrco_Leaf_FDR[i,]$C_delta_Parents - Scrco_Leaf_FDR[i,]$C_delta_CASI_co) / rowMaxs(as.matrix(Scrco_Leaf_FDR[i, c("C_delta_Parents","C_delta_CASI_co")]))
}

head(Scrco_Leaf_FDR)

# Root
head(Scrco_Root_FDR)

# Within dataframe of differently expressed between parental species. CR CO
# Absolute differences between parents in expression.
Scrco_Root_FDR$C_delta_Parents <- abs(Scrco_Root_FDR$mean_CR_R - Scrco_Root_FDR$mean_CO_R)

# Absolute differences between Subgenomes in expression.
Scrco_Root_FDR$C_delta_Cbp <- abs(Scrco_Root_FDR$mean_Cbp_Cg_R - Scrco_Root_FDR$mean_Cbp_Co_R) # Cbp all
Scrco_Root_FDR$C_delta_ASI <- abs(Scrco_Root_FDR$mean_ASI_Cg_R - Scrco_Root_FDR$mean_ASI_Co_R) # Cbp-ASI
Scrco_Root_FDR$C_delta_EUR <- abs(Scrco_Root_FDR$mean_EUR_Cg_R - Scrco_Root_FDR$mean_EUR_Co_R) # Cbp-EUR
Scrco_Root_FDR$C_delta_ME <- abs(Scrco_Root_FDR$mean_ME_Cg_R - Scrco_Root_FDR$mean_ME_Co_R) # Cbp-ME
Scrco_Root_FDR$C_delta_CASI <- abs(Scrco_Root_FDR$mean_CASI_Cg_R - Scrco_Root_FDR$mean_CASI_Co_R) # Cbp-CASI

# each subgenome and the opposite parental species
Scrco_Root_FDR$C_delta_Cbp_cg <- abs(Scrco_Root_FDR$mean_Cbp_Cg_R - Scrco_Root_FDR$mean_CO_R) # Cbp all
Scrco_Root_FDR$C_delta_Cbp_co <- abs(Scrco_Root_FDR$mean_Cbp_Co_R - Scrco_Root_FDR$mean_CR_R)

Scrco_Root_FDR$C_delta_ASI_cg <- abs(Scrco_Root_FDR$mean_ASI_Cg_R - Scrco_Root_FDR$mean_CO_R) # ASI
Scrco_Root_FDR$C_delta_ASI_co <- abs(Scrco_Root_FDR$mean_ASI_Co_R - Scrco_Root_FDR$mean_CR_R)

Scrco_Root_FDR$C_delta_EUR_cg <- abs(Scrco_Root_FDR$mean_EUR_Cg_R - Scrco_Root_FDR$mean_CO_R) # EUR
Scrco_Root_FDR$C_delta_EUR_co <- abs(Scrco_Root_FDR$mean_EUR_Co_R - Scrco_Root_FDR$mean_CR_R)

Scrco_Root_FDR$C_delta_ME_cg <- abs(Scrco_Root_FDR$mean_ME_Cg_R - Scrco_Root_FDR$mean_CO_R) # ME
Scrco_Root_FDR$C_delta_ME_co <- abs(Scrco_Root_FDR$mean_ME_Co_R - Scrco_Root_FDR$mean_CR_R)

Scrco_Root_FDR$C_delta_CASI_cg <- abs(Scrco_Root_FDR$mean_CASI_Cg_R - Scrco_Root_FDR$mean_CO_R) # CASI
Scrco_Root_FDR$C_delta_CASI_co <- abs(Scrco_Root_FDR$mean_CASI_Co_R - Scrco_Root_FDR$mean_CR_R)

head(Scrco_Root_FDR)

# define columns with empty value.
Scrco_Root_FDR$C.Cbp <- NA
Scrco_Root_FDR$C.Cbp_cg <- NA
Scrco_Root_FDR$C.Cbp_co <- NA

Scrco_Root_FDR$C.ASI <- NA
Scrco_Root_FDR$C.ASI_cg <- NA
Scrco_Root_FDR$C.ASI_co <- NA

Scrco_Root_FDR$C.EUR <- NA
Scrco_Root_FDR$C.EUR_cg<- NA
Scrco_Root_FDR$C.EUR_co <- NA

Scrco_Root_FDR$C.ME <- NA
Scrco_Root_FDR$C.ME_cg <- NA
Scrco_Root_FDR$C.ME_co <- NA

Scrco_Root_FDR$C.CASI <- NA
Scrco_Root_FDR$C.CASI_cg <- NA
Scrco_Root_FDR$C.CASI_co <- NA

# Perform the C index calculation
for(i in c(1:dim(Scrco_Root_FDR)[1])){
  Scrco_Root_FDR[i,]$C.Cbp <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_Cbp) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_Cbp")]))
  Scrco_Root_FDR[i,]$C.Cbp_cg <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_Cbp_cg) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_Cbp_cg")]))
  Scrco_Root_FDR[i,]$C.Cbp_co <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_Cbp_co) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_Cbp_co")]))
  
  Scrco_Root_FDR[i,]$C.ASI <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ASI) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ASI")]))
  Scrco_Root_FDR[i,]$C.ASI_cg <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ASI_cg) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ASI_cg")]))
  Scrco_Root_FDR[i,]$C.ASI_co <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ASI_co) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ASI_co")]))
  
  Scrco_Root_FDR[i,]$C.EUR <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_EUR) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_EUR")]))
  Scrco_Root_FDR[i,]$C.EUR_cg <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_EUR_cg) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_EUR_cg")]))
  Scrco_Root_FDR[i,]$C.EUR_co <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_EUR_co) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_EUR_co")]))
  
  Scrco_Root_FDR[i,]$C.ME <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ME) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ME")]))
  Scrco_Root_FDR[i,]$C.ME_cg <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ME_cg) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ME_cg")]))
  Scrco_Root_FDR[i,]$C.ME_co <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_ME_co) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_ME_co")]))
  
  Scrco_Root_FDR[i,]$C.CASI <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_CASI) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_CASI")]))
  Scrco_Root_FDR[i,]$C.CASI_cg <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_CASI_cg) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_CASI_cg")]))
  Scrco_Root_FDR[i,]$C.CASI_co <- (Scrco_Root_FDR[i,]$C_delta_Parents - Scrco_Root_FDR[i,]$C_delta_CASI_co) / rowMaxs(as.matrix(Scrco_Root_FDR[i, c("C_delta_Parents","C_delta_CASI_co")]))
}

head(Scrco_Root_FDR)

# Output data with values of S and C indices.
# write.table(Scrco_Flower_FDR, file = "OutputData/S and C indices in Flower FDR genes.txt", quote = F, sep = "\t", row.names = F)
# write.table(Scrco_Leaf_FDR, file = "OutputData/S and C indices in Leaf FDR genes.txt", quote = F, sep = "\t", row.names = F)
# write.table(Scrco_Root_FDR, file = "OutputData/S and C indices in Root FDR genes.txt", quote = F, sep = "\t", row.names = F)

# Create a new dataframe of S for plot
S_index <- c(rep(c("delta_S","S_Co", "S_Cg"), 15)) # 5 Cbp subgroups (all + asi + eur + me + casi), 3 tissues (F,L,R)
S_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 3)
S_tissue <- c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15))
# Count median S value for each subpopulation. Median value
S_value1 <- c(median(Scrco_Flower_FDR$delta_Scrco_Cbp), # Flower all  ##### Flower
           median(Scrco_Flower_FDR$Scrco_Cbp_Co),
           median(Scrco_Flower_FDR$Scrco_Cbp_Cg),
           median(Scrco_Flower_FDR$delta_Scrco_ASI), # Flower ASI
           median(Scrco_Flower_FDR$Scrco_ASI_Co),
           median(Scrco_Flower_FDR$Scrco_ASI_Cg),
           median(Scrco_Flower_FDR$delta_Scrco_EUR), # Flower EUR
           median(Scrco_Flower_FDR$Scrco_EUR_Co),
           median(Scrco_Flower_FDR$Scrco_EUR_Cg),
           median(Scrco_Flower_FDR$delta_Scrco_ME), # Flower ME
           median(Scrco_Flower_FDR$Scrco_ME_Co),
           median(Scrco_Flower_FDR$Scrco_ME_Cg),
           median(Scrco_Flower_FDR$delta_Scrco_CASI), # Flower CASI
           median(Scrco_Flower_FDR$Scrco_CASI_Co),
           median(Scrco_Flower_FDR$Scrco_CASI_Cg),
           median(Scrco_Leaf_FDR$delta_Scrco_Cbp), # Leaf all  ##### Leaf
           median(Scrco_Leaf_FDR$Scrco_Cbp_Co),
           median(Scrco_Leaf_FDR$Scrco_Cbp_Cg),
           median(Scrco_Leaf_FDR$delta_Scrco_ASI), # Leaf ASI
           median(Scrco_Leaf_FDR$Scrco_ASI_Co),
           median(Scrco_Leaf_FDR$Scrco_ASI_Cg),
           median(Scrco_Leaf_FDR$delta_Scrco_EUR), # Leaf EUR
           median(Scrco_Leaf_FDR$Scrco_EUR_Co),
           median(Scrco_Leaf_FDR$Scrco_EUR_Cg),
           median(Scrco_Leaf_FDR$delta_Scrco_ME), # Leaf ME
           median(Scrco_Leaf_FDR$Scrco_ME_Co),
           median(Scrco_Leaf_FDR$Scrco_ME_Cg),
           median(Scrco_Leaf_FDR$delta_Scrco_CASI), # Leaf CASI
           median(Scrco_Leaf_FDR$Scrco_CASI_Co),
           median(Scrco_Leaf_FDR$Scrco_CASI_Cg),
           median(Scrco_Root_FDR$delta_Scrco_Cbp), # Root all  ##### Root
           median(Scrco_Root_FDR$Scrco_Cbp_Co),
           median(Scrco_Root_FDR$Scrco_Cbp_Cg),
           median(Scrco_Root_FDR$delta_Scrco_ASI), # Root ASI
           median(Scrco_Root_FDR$Scrco_ASI_Co),
           median(Scrco_Root_FDR$Scrco_ASI_Cg),
           median(Scrco_Root_FDR$delta_Scrco_EUR), # Root EUR
           median(Scrco_Root_FDR$Scrco_EUR_Co),
           median(Scrco_Root_FDR$Scrco_EUR_Cg),
           median(Scrco_Root_FDR$delta_Scrco_ME), # Root ME
           median(Scrco_Root_FDR$Scrco_ME_Co),
           median(Scrco_Root_FDR$Scrco_ME_Cg),
           median(Scrco_Root_FDR$delta_Scrco_CASI), # Root CASI
           median(Scrco_Root_FDR$Scrco_CASI_Co),
           median(Scrco_Root_FDR$Scrco_CASI_Cg)
           )
S_value1

S_value2 <- c(abs(median(Scrco_Flower_FDR$Scrco_Cbp_Co)) - abs(median(Scrco_Flower_FDR$Scrco_Cbp_Cg)), # Flower all  ##### Flower
              median(Scrco_Flower_FDR$Scrco_Cbp_Co),
              median(Scrco_Flower_FDR$Scrco_Cbp_Cg),
              abs(median(Scrco_Flower_FDR$Scrco_ASI_Co)) - abs(median(Scrco_Flower_FDR$Scrco_ASI_Cg)), # Flower ASI
              median(Scrco_Flower_FDR$Scrco_ASI_Co),
              median(Scrco_Flower_FDR$Scrco_ASI_Cg),
              abs(median(Scrco_Flower_FDR$Scrco_EUR_Co)) - abs(median(Scrco_Flower_FDR$Scrco_EUR_Cg)), # Flower EUR
              median(Scrco_Flower_FDR$Scrco_EUR_Co),
              median(Scrco_Flower_FDR$Scrco_EUR_Cg),
              abs(median(Scrco_Flower_FDR$Scrco_ME_Co)) - abs(median(Scrco_Flower_FDR$Scrco_ME_Cg)), # Flower ME
              median(Scrco_Flower_FDR$Scrco_ME_Co),
              median(Scrco_Flower_FDR$Scrco_ME_Cg),
              abs(median(Scrco_Flower_FDR$Scrco_CASI_Co)) - abs(median(Scrco_Flower_FDR$Scrco_CASI_Cg)), # Flower CASI
              median(Scrco_Flower_FDR$Scrco_CASI_Co),
              median(Scrco_Flower_FDR$Scrco_CASI_Cg),
              abs(median(Scrco_Leaf_FDR$Scrco_Cbp_Co)) - abs(median(Scrco_Leaf_FDR$Scrco_Cbp_Cg)), # Leaf all  ##### Leaf
              median(Scrco_Leaf_FDR$Scrco_Cbp_Co),
              median(Scrco_Leaf_FDR$Scrco_Cbp_Cg),
              abs(median(Scrco_Leaf_FDR$Scrco_ASI_Co)) - abs(median(Scrco_Leaf_FDR$Scrco_ASI_Cg)), # Leaf ASI
              median(Scrco_Leaf_FDR$Scrco_ASI_Co),
              median(Scrco_Leaf_FDR$Scrco_ASI_Cg),
              abs(median(Scrco_Leaf_FDR$Scrco_EUR_Co)) - abs(median(Scrco_Leaf_FDR$Scrco_EUR_Cg)), # Leaf EUR
              median(Scrco_Leaf_FDR$Scrco_EUR_Co),
              median(Scrco_Leaf_FDR$Scrco_EUR_Cg),
              abs(median(Scrco_Leaf_FDR$Scrco_ME_Co)) - abs(median(Scrco_Leaf_FDR$Scrco_ME_Cg)), # Leaf ME
              median(Scrco_Leaf_FDR$Scrco_ME_Co),
              median(Scrco_Leaf_FDR$Scrco_ME_Cg),
              abs(median(Scrco_Leaf_FDR$Scrco_CASI_Co)) - abs(median(Scrco_Leaf_FDR$Scrco_CASI_Cg)), # Leaf CASI
              median(Scrco_Leaf_FDR$Scrco_CASI_Co),
              median(Scrco_Leaf_FDR$Scrco_CASI_Cg),
              abs(median(Scrco_Root_FDR$Scrco_Cbp_Co)) - abs(median(Scrco_Root_FDR$Scrco_Cbp_Cg)), # Root all  ##### Root
              median(Scrco_Root_FDR$Scrco_Cbp_Co),
              median(Scrco_Root_FDR$Scrco_Cbp_Cg),
              abs(median(Scrco_Root_FDR$Scrco_ASI_Co)) - abs(median(Scrco_Root_FDR$Scrco_ASI_Cg)), # Root ASI
              median(Scrco_Root_FDR$Scrco_ASI_Co),
              median(Scrco_Root_FDR$Scrco_ASI_Cg),
              abs(median(Scrco_Root_FDR$Scrco_EUR_Co)) - abs(median(Scrco_Root_FDR$Scrco_EUR_Cg)), # Root EUR
              median(Scrco_Root_FDR$Scrco_EUR_Co),
              median(Scrco_Root_FDR$Scrco_EUR_Cg),
              abs(median(Scrco_Root_FDR$Scrco_ME_Co)) - abs(median(Scrco_Root_FDR$Scrco_ME_Cg)), # Root ME
              median(Scrco_Root_FDR$Scrco_ME_Co),
              median(Scrco_Root_FDR$Scrco_ME_Cg),
              abs(median(Scrco_Root_FDR$Scrco_CASI_Co)) - abs(median(Scrco_Root_FDR$Scrco_CASI_Cg)), # Root CASI
              median(Scrco_Root_FDR$Scrco_CASI_Co),
              median(Scrco_Root_FDR$Scrco_CASI_Cg)
)
S_value2

# Combine all vectors to a dataframe
S_plot <- cbind(S_index, S_population, S_tissue, S_value2)
# change the format of the dataframe
S_plot <- print.data.frame(data.frame(S_plot), quote=FALSE)
str(S_plot)
S_plot$S_index <- factor(S_plot$S_index, levels = c("S_Co","S_Cg","delta_S"))
S_plot$S_population <- factor(S_plot$S_population, levels = c("Cbp","Cbp_ASI","Cbp_EUR","Cbp_ME","Cbp_CASI"))
S_plot$S_tissue <- factor(S_plot$S_tissue, levels = c("Flower","Leaf","Root"))
S_plot$S_value <- as.numeric(S_plot$S_value2)
str(S_plot)

# Plot
library(ggplot2)
p1 <- ggplot(S_plot, aes(x=S_tissue, y=S_value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=S_index, shape=S_index),size=2) +
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  #scale_y_continuous(limits=c(-0.15, 0.15)) +
  facet_grid(. ~ S_population) + 
  labs(x ="Tissues", y = "Median Si") +
  theme_bw()

ggplot(S_plot, aes(x=S_tissue, y=S_value, color=S_index, shape=S_index)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(size=3) +
  #scale_y_continuous(limits=c(-0.15, 0.15)) +
  facet_grid(. ~ S_population) + 
  labs(x ="Tissues", y = "Median Si") +
  theme_bw()
  #theme(panel.background = element_rect(fill = "white", colour = "grey50"))

library(ggpubr)
ggscatter(S_plot, x = "S_tissue", y = "S_value", color = "S_index", shape = "S_index",
          ylim = c(-0.3,0.3))

# Create a new dataframe of C for plot
C_index <- c(rep(c("Cbp","Cbp_co", "Cbp_cg"), 15)) # 5 Cbp subgroups (all + asi + eur + me + casi), 3 tissues (F,L,R)
C_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 3)
C_tissue <- c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15))
# Count proportion of C > 0
C_value <- c(nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.Cbp > 0, ]) / nrow(Scrco_Flower_FDR), # Flower all  ##### Flower
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.Cbp_co > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.Cbp_cg > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ASI > 0, ]) / nrow(Scrco_Flower_FDR), # Flower ASI
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ASI_co > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ASI_cg > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.EUR > 0, ]) / nrow(Scrco_Flower_FDR), # Flower EUR
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.EUR_co > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.EUR_cg > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ME > 0, ]) / nrow(Scrco_Flower_FDR), # Flower ME
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ME_co > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.ME_cg > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.CASI > 0, ]) / nrow(Scrco_Flower_FDR), # Flower CASI
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.CASI_co > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Flower_FDR[Scrco_Flower_FDR$C.CASI_cg > 0, ]) / nrow(Scrco_Flower_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.Cbp > 0, ]) / nrow(Scrco_Leaf_FDR), # Leaf all  ##### Leaf
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.Cbp_co > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.Cbp_cg > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ASI > 0, ]) / nrow(Scrco_Leaf_FDR), # Leaf ASI
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ASI_co > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ASI_cg > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.EUR > 0, ]) / nrow(Scrco_Leaf_FDR), # Leaf EUR
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.EUR_co > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.EUR_cg > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ME > 0, ]) / nrow(Scrco_Leaf_FDR), # Leaf ME
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ME_co > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.ME_cg > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.CASI > 0, ]) / nrow(Scrco_Leaf_FDR), # Leaf CASI
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.CASI_co > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Leaf_FDR[Scrco_Leaf_FDR$C.CASI_cg > 0, ]) / nrow(Scrco_Leaf_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.Cbp > 0, ]) / nrow(Scrco_Root_FDR), # Root all  ##### Root
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.Cbp_co > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.Cbp_cg > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ASI > 0, ]) / nrow(Scrco_Root_FDR), # Root ASI
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ASI_co > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ASI_cg > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.EUR > 0, ]) / nrow(Scrco_Root_FDR), # Root EUR
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.EUR_co > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.EUR_cg > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ME > 0, ]) / nrow(Scrco_Root_FDR), # Root ME
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ME_co > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.ME_cg > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.CASI > 0, ]) / nrow(Scrco_Root_FDR), # Root CASI
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.CASI_co > 0, ]) / nrow(Scrco_Root_FDR),
           nrow(Scrco_Root_FDR[Scrco_Root_FDR$C.CASI_cg > 0, ]) / nrow(Scrco_Root_FDR)
)
C_value

C_plot <- cbind(C_index, C_population, C_tissue, C_value)
# change the format of the dataframe
C_plot <- print.data.frame(data.frame(C_plot), quote=FALSE)
str(C_plot)
C_plot$C_index <- factor(C_plot$C_index, levels = c("Cbp","Cbp_co","Cbp_cg"))
C_plot$C_population <- factor(C_plot$C_population, levels = c("Cbp","Cbp_ASI","Cbp_EUR","Cbp_ME","Cbp_CASI"))
C_plot$C_tissue <- factor(C_plot$C_tissue, levels = c("Flower","Leaf","Root"))
C_plot$C_value <- as.numeric(C_plot$C_value)
str(C_plot)

# Plot
library(ggplot2)

p2 <- ggplot(C_plot, aes(x = C_tissue, y = C_value, colour= C_index)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey31") +
  geom_point(aes(shape=C_index), position=position_dodge(width=0.5), size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(3, 4, 8))+
  scale_color_manual(values=c("black","dodgerblue", "brown1"))+
  #scale_color_manual(values=c('#E69F00','#56B4E9',"tomato" ))+
  #scale_y_continuous(limits=c(0.4, 1.0), breaks = seq(0.4, 1.0, by = 0.1)) +
  facet_grid(. ~ C_population) + 
  labs(x ="Tissues", y = "Proportion of Ci > 0") +
  theme_bw()

library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(2,1))

###### Distribution of S index across tissues. #########
# flower
head(Scrco_Flower_FDR)

p1 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_Cbp_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_cg", y = "Density", title = "Flower All") +
  theme_classic()

p2 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_Cbp_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_co", y = "Density", title = "") +
  theme_classic()

p3 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_ASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_cg", y = "Density", title = "Flower ASI") +
  theme_bw()

p4 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_ASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_co", y = "Density", title = "") +
  theme_bw()

p5 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_EUR_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_cg", y = "Density", title = "Flower EUR") +
  theme_bw()

p6 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_EUR_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_co", y = "Density", title = "Flower EUR") +
  theme_bw()

p7 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_ME_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_cg", y = "Density", title = "Flower ME") +
  theme_bw()

p8 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_ME_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_co", y = "Density", title = "") +
  theme_bw()

p9 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_CASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_cg", y = "Density", title = "Flower CASI") +
  theme_bw()

p10 <- ggplot(Scrco_Flower_FDR, aes(x=Scrco_CASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_co", y = "Density", title = "") +
  theme_bw()

p11 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_Cbp_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_cg", y = "Density", title = "Leaf All") +
  theme_classic()

p12 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_Cbp_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_co", y = "Density", title = "") +
  theme_classic()

p13 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_ASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_cg", y = "Density", title = "Leaf ASI") +
  theme_bw()

p14 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_ASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_co", y = "Density", title = "") +
  theme_bw()

p15 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_EUR_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_cg", y = "Density", title = "Leaf EUR") +
  theme_bw()

p16 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_EUR_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_co", y = "Density", title = "Leaf EUR") +
  theme_bw()

p17 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_ME_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_cg", y = "Density", title = "Leaf ME") +
  theme_bw()

p18 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_ME_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_co", y = "Density", title = "") +
  theme_bw()

p19 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_CASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_cg", y = "Density", title = "Leaf CASI") +
  theme_bw()

p20 <- ggplot(Scrco_Leaf_FDR, aes(x=Scrco_CASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_co", y = "Density", title = "") +
  theme_bw()

p21 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_Cbp_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_cg", y = "Density", title = "Root All") +
  theme_classic()

p22 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_Cbp_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_Cbp_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S Cbp_co", y = "Density", title = "") +
  theme_classic()

p23 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_ASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_cg", y = "Density", title = "Root ASI") +
  theme_bw()

p24 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_ASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ASI_co", y = "Density", title = "") +
  theme_bw()

p25 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_EUR_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_cg", y = "Density", title = "Root EUR") +
  theme_bw()

p26 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_EUR_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_EUR_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S EUR_co", y = "Density", title = "Root EUR") +
  theme_bw()

p27 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_ME_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_cg", y = "Density", title = "Root ME") +
  theme_bw()

p28 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_ME_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_ME_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S ME_co", y = "Density", title = "") +
  theme_bw()

p29 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_CASI_Cg)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Cg)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_cg", y = "Density", title = "Root CASI") +
  theme_bw()

p30 <- ggplot(Scrco_Root_FDR, aes(x=Scrco_CASI_Co)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 50) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(Scrco_CASI_Co)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CASI_co", y = "Density", title = "") +
  theme_bw()


library("grid")
grid.newpage()
pushViewport(viewport(layout = grid.layout(6,5)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p3, vp = vplayout(1,3))
print(p4, vp = vplayout(1,4))
print(p5, vp = vplayout(1,5))
print(p11, vp = vplayout(2,1))
print(p12, vp = vplayout(2,2))
print(p13, vp = vplayout(2,3))
print(p14, vp = vplayout(2,4))
print(p15, vp = vplayout(2,5))
print(p21, vp = vplayout(3,1))
print(p22, vp = vplayout(3,2))
print(p23, vp = vplayout(3,3))
print(p24, vp = vplayout(3,4))
print(p25, vp = vplayout(3,5))
print(p6, vp = vplayout(4,1))
print(p7, vp = vplayout(4,2))
print(p8, vp = vplayout(4,3))
print(p9, vp = vplayout(4,4))
print(p10, vp = vplayout(4,5))
print(p16, vp = vplayout(5,1))
print(p17, vp = vplayout(5,2))
print(p18, vp = vplayout(5,3))
print(p19, vp = vplayout(5,4))
print(p20, vp = vplayout(5,5))
print(p26, vp = vplayout(6,1))
print(p27, vp = vplayout(6,2))
print(p28, vp = vplayout(6,3))
print(p29, vp = vplayout(6,4))
print(p30, vp = vplayout(6,5))


grid.newpage()
pushViewport(viewport(layout = grid.layout(3,2)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p11, vp = vplayout(2,1))
print(p12, vp = vplayout(2,2))
print(p21, vp = vplayout(3,1))
print(p22, vp = vplayout(3,2))




