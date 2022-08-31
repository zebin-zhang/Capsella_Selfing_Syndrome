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

CPM <- read.table("InputData/Capsella_ReadsCount_CPM_Filtered.txt", header = T, sep = "\t")
head(CPM)
dim(CPM)
CPM[CPM == 0] <- NA # Substitute 0 to NA.

logCPM <- read.table("InputData/Capsella_ReadsCount_logCPM_Filtered.txt", header = T, sep = "\t")
head(logCPM)
dim(logCPM)
CPM[logCPM == 0] <- NA 

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
######################  CPM   ########################
# Similarity index
# CPM
# Create a new dataframe for Flower
CPM_S_Flower <- data.frame(Gene=rownames(CPM))
CPM_S_Flower$mean_CG_F <- rowMeans(CPM[,CG_F], na.rm = T) # Average value of CG expression in Flower
CPM_S_Flower$mean_CR_F <- rowMeans(CPM[,CR_F], na.rm = T)
CPM_S_Flower$mean_CO_F <- rowMeans(CPM[,CO_F], na.rm = T)

CPM_S_Flower$mean_Cbp_Cg_F <- rowMeans(CPM[,Cbp_Cg_F], na.rm = T)
CPM_S_Flower$mean_Cbp_Co_F <- rowMeans(CPM[,Cbp_Co_F], na.rm = T)

CPM_S_Flower$mean_ASI_Cg_F <- rowMeans(CPM[,ASI_Cg_F], na.rm = T)
CPM_S_Flower$mean_EUR_Cg_F <- rowMeans(CPM[,EUR_Cg_F], na.rm = T)
CPM_S_Flower$mean_ME_Cg_F <- rowMeans(CPM[,ME_Cg_F], na.rm = T)
CPM_S_Flower$mean_CASI_Cg_F <- rowMeans(CPM[,CASI_Cg_F], na.rm = T)

CPM_S_Flower$mean_ASI_Co_F <- rowMeans(CPM[,ASI_Co_F], na.rm = T)
CPM_S_Flower$mean_EUR_Co_F <- rowMeans(CPM[,EUR_Co_F], na.rm = T)
CPM_S_Flower$mean_ME_Co_F <- rowMeans(CPM[,ME_Co_F], na.rm = T)
CPM_S_Flower$mean_CASI_Co_F <- rowMeans(CPM[,CASI_Co_F], na.rm = T)

head(CPM[,CG_F])
head(CPM[,EUR_Cg_F])
head(CPM_S_Flower)

# Average expression level in parents of CR & CO.
CPM_S_Flower$Ave_CRCO <- rowMeans(CPM_S_Flower[,c("mean_CR_F","mean_CO_F")])
CPM_S_Flower$Ave_CGCO <- rowMeans(CPM_S_Flower[,c("mean_CG_F","mean_CO_F")])
CPM_S_Flower$Ave_CGCR <- rowMeans(CPM_S_Flower[,c("mean_CG_F","mean_CR_F")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
CPM_S_Flower$CPM_Scrco_Cbp_Cg <- (CPM_S_Flower$mean_Cbp_Cg_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO # Subgenome Cg 

CPM_S_Flower$CPM_Scrco_ASI_Cg <- (CPM_S_Flower$mean_ASI_Cg_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_EUR_Cg <- (CPM_S_Flower$mean_EUR_Cg_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_ME_Cg <- (CPM_S_Flower$mean_ME_Cg_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_CASI_Cg <- (CPM_S_Flower$mean_CASI_Cg_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO

CPM_S_Flower$CPM_Scrco_Cbp_Co <- (CPM_S_Flower$mean_Cbp_Co_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO # Subgenome Co

CPM_S_Flower$CPM_Scrco_ASI_Co <- (CPM_S_Flower$mean_ASI_Co_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_EUR_Co <- (CPM_S_Flower$mean_EUR_Co_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_ME_Co <- (CPM_S_Flower$mean_ME_Co_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO
CPM_S_Flower$CPM_Scrco_CASI_Co <- (CPM_S_Flower$mean_CASI_Co_F - CPM_S_Flower$Ave_CRCO) / CPM_S_Flower$Ave_CRCO

head(CPM_S_Flower)
# in CG vs CO
CPM_S_Flower$CPM_Scgco_Cbp_Cg <- (CPM_S_Flower$mean_Cbp_Cg_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO  # Subgenome Cg 

CPM_S_Flower$CPM_Scgco_ASI_Cg <- (CPM_S_Flower$mean_ASI_Cg_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_EUR_Cg <- (CPM_S_Flower$mean_EUR_Cg_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_ME_Cg <- (CPM_S_Flower$mean_ME_Cg_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_CASI_Cg <- (CPM_S_Flower$mean_CASI_Cg_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO

CPM_S_Flower$CPM_Scgco_Cbp_Co <- (CPM_S_Flower$mean_Cbp_Co_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO  # Subgenome Co

CPM_S_Flower$CPM_Scgco_ASI_Co <- (CPM_S_Flower$mean_ASI_Co_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_EUR_Co <- (CPM_S_Flower$mean_EUR_Co_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_ME_Co <- (CPM_S_Flower$mean_ME_Co_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO
CPM_S_Flower$CPM_Scgco_CASI_Co <- (CPM_S_Flower$mean_CASI_Co_F - CPM_S_Flower$Ave_CGCO) / CPM_S_Flower$Ave_CGCO

head(CPM_S_Flower)
# in CG vs CR
CPM_S_Flower$CPM_Scgcr_Cbp_Cg <- (CPM_S_Flower$mean_Cbp_Cg_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR  # Subgenome Cg 

CPM_S_Flower$CPM_Scgcr_ASI_Cg <- (CPM_S_Flower$mean_ASI_Cg_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_EUR_Cg <- (CPM_S_Flower$mean_EUR_Cg_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_ME_Cg <- (CPM_S_Flower$mean_ME_Cg_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_CASI_Cg <- (CPM_S_Flower$mean_CASI_Cg_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR

CPM_S_Flower$CPM_Scgcr_Cbp_Co <- (CPM_S_Flower$mean_Cbp_Co_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR  # Subgenome Co

CPM_S_Flower$CPM_Scgcr_ASI_Co <- (CPM_S_Flower$mean_ASI_Co_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_EUR_Co <- (CPM_S_Flower$mean_EUR_Co_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_ME_Co <- (CPM_S_Flower$mean_ME_Co_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR
CPM_S_Flower$CPM_Scgcr_CASI_Co <- (CPM_S_Flower$mean_CASI_Co_F - CPM_S_Flower$Ave_CGCR) / CPM_S_Flower$Ave_CGCR

head(CPM_S_Flower)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_Cbp_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_Cbp_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_EUR_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_EUR_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ME_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ME_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_CASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_CASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_Cbp_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_Cbp_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ASI_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_EUR_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_EUR_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ME_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_ME_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_CASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CR_F > CPM_S_Flower$mean_CO_F),]$CPM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_Cbp_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_Cbp_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_EUR_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_EUR_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ME_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ME_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_CASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_CASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_Cbp_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_Cbp_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ASI_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_EUR_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_EUR_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ME_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_ME_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_CASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CO_F),]$CPM_Scgco_CASI_Co * (-1)



# CG vs CR, CG > CR.

CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_Cbp_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_Cbp_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_EUR_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_EUR_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ME_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ME_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_CASI_Cg <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_CASI_Cg * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_Cbp_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_Cbp_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ASI_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_EUR_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_EUR_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ME_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_ME_Co * (-1)
CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_CASI_Co <- CPM_S_Flower[(CPM_S_Flower$mean_CG_F > CPM_S_Flower$mean_CR_F),]$CPM_Scgcr_CASI_Co * (-1)

# Calculate delta S by comparing the absolute values of CPM_S_Cbp_Cr and CPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
CPM_S_Flower$delta_CPM_Scrco_Cbp <- abs(CPM_S_Flower$CPM_Scrco_Cbp_Cg) - abs(CPM_S_Flower$CPM_Scrco_Cbp_Co)
CPM_S_Flower$delta_CPM_Scrco_ASI <- abs(CPM_S_Flower$CPM_Scrco_ASI_Cg) - abs(CPM_S_Flower$CPM_Scrco_ASI_Co)
CPM_S_Flower$delta_CPM_Scrco_EUR <- abs(CPM_S_Flower$CPM_Scrco_EUR_Cg) - abs(CPM_S_Flower$CPM_Scrco_EUR_Co)
CPM_S_Flower$delta_CPM_Scrco_ME <- abs(CPM_S_Flower$CPM_Scrco_ME_Cg) - abs(CPM_S_Flower$CPM_Scrco_ME_Co)
CPM_S_Flower$delta_CPM_Scrco_CASI <- abs(CPM_S_Flower$CPM_Scrco_CASI_Cg) - abs(CPM_S_Flower$CPM_Scrco_CASI_Co)
head(CPM_S_Flower)

# CG vs CO
CPM_S_Flower$delta_CPM_Scgco_Cbp <- abs(CPM_S_Flower$CPM_Scgco_Cbp_Cg) - abs(CPM_S_Flower$CPM_Scgco_Cbp_Co)
CPM_S_Flower$delta_CPM_Scgco_ASI <- abs(CPM_S_Flower$CPM_Scgco_ASI_Cg) - abs(CPM_S_Flower$CPM_Scgco_ASI_Co)
CPM_S_Flower$delta_CPM_Scgco_EUR <- abs(CPM_S_Flower$CPM_Scgco_EUR_Cg) - abs(CPM_S_Flower$CPM_Scgco_EUR_Co)
CPM_S_Flower$delta_CPM_Scgco_ME <- abs(CPM_S_Flower$CPM_Scgco_ME_Cg) - abs(CPM_S_Flower$CPM_Scgco_ME_Co)
CPM_S_Flower$delta_CPM_Scgco_CASI <- abs(CPM_S_Flower$CPM_Scgco_CASI_Cg) - abs(CPM_S_Flower$CPM_Scgco_CASI_Co)
head(CPM_S_Flower)

# CG vs CR
CPM_S_Flower$delta_CPM_Scgcr_Cbp <- abs(CPM_S_Flower$CPM_Scgcr_Cbp_Cg) - abs(CPM_S_Flower$CPM_Scgcr_Cbp_Co)
CPM_S_Flower$delta_CPM_Scgcr_ASI <- abs(CPM_S_Flower$CPM_Scgcr_ASI_Cg) - abs(CPM_S_Flower$CPM_Scgcr_ASI_Co)
CPM_S_Flower$delta_CPM_Scgcr_EUR <- abs(CPM_S_Flower$CPM_Scgcr_EUR_Cg) - abs(CPM_S_Flower$CPM_Scgcr_EUR_Co)
CPM_S_Flower$delta_CPM_Scgcr_ME <- abs(CPM_S_Flower$CPM_Scgcr_ME_Cg) - abs(CPM_S_Flower$CPM_Scgcr_ME_Co)
CPM_S_Flower$delta_CPM_Scgcr_CASI <- abs(CPM_S_Flower$CPM_Scgcr_CASI_Cg) - abs(CPM_S_Flower$CPM_Scgcr_CASI_Co)
head(CPM_S_Flower)
# load genes significantly expressed between Cr Co in flower.
DE.F.crco <- read.table("OutputData/DE_CrCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.crco)
dim(DE.F.crco)

DE.F.cgco <- read.table("OutputData/DE_CgCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgco)
dim(DE.F.cgco)

DE.F.cgcr <- read.table("OutputData/DE_CgCr_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgcr)
dim(DE.F.cgcr)

# keep genes with FDR ≤ 0.01
DE.F.crco.keep <- CPM_S_Flower$Gene %in% row.names(DE.F.crco)
CPM_S_Flower$Gene %in% row.names(DE.F.crco)
DE.F.cgco.keep <- CPM_S_Flower$Gene %in% row.names(DE.F.cgco)
CPM_S_Flower$Gene %in% row.names(DE.F.cgco)
DE.F.cgcr.keep <- CPM_S_Flower$Gene %in% row.names(DE.F.cgcr)
CPM_S_Flower$Gene %in% row.names(DE.F.cgcr)

CPM_Scrco_Flower_FDR <- CPM_S_Flower[DE.F.crco.keep,]
head(CPM_Scrco_Flower_FDR)
dim(CPM_Scrco_Flower_FDR)
CPM_Scgco_Flower_FDR <- CPM_S_Flower[DE.F.cgco.keep,]
head(CPM_Scgco_Flower_FDR)
dim(CPM_Scgco_Flower_FDR)
CPM_Scgcr_Flower_FDR <- CPM_S_Flower[DE.F.cgcr.keep,]
head(CPM_Scgcr_Flower_FDR)
dim(CPM_Scgcr_Flower_FDR)

##################################

# Create a new dataframe for Leaf
CPM_S_Leaf <- data.frame(Gene=rownames(CPM))
CPM_S_Leaf$mean_CG_L <- rowMeans(CPM[,CG_L], na.rm = T) # Average value of CG expression in Leaf
CPM_S_Leaf$mean_CR_L <- rowMeans(CPM[,CR_L], na.rm = T)
CPM_S_Leaf$mean_CO_L <- rowMeans(CPM[,CO_L], na.rm = T)

CPM_S_Leaf$mean_Cbp_Cg_L <- rowMeans(CPM[,Cbp_Cg_L], na.rm = T)
CPM_S_Leaf$mean_Cbp_Co_L <- rowMeans(CPM[,Cbp_Co_L], na.rm = T)

CPM_S_Leaf$mean_ASI_Cg_L <- rowMeans(CPM[,ASI_Cg_L], na.rm = T)
CPM_S_Leaf$mean_EUR_Cg_L <- rowMeans(CPM[,EUR_Cg_L], na.rm = T)
CPM_S_Leaf$mean_ME_Cg_L <- rowMeans(CPM[,ME_Cg_L], na.rm = T)
CPM_S_Leaf$mean_CASI_Cg_L <- rowMeans(CPM[,CASI_Cg_L], na.rm = T)

CPM_S_Leaf$mean_ASI_Co_L <- rowMeans(CPM[,ASI_Co_L], na.rm = T)
CPM_S_Leaf$mean_EUR_Co_L <- rowMeans(CPM[,EUR_Co_L], na.rm = T)
CPM_S_Leaf$mean_ME_Co_L <- rowMeans(CPM[,ME_Co_L], na.rm = T)
CPM_S_Leaf$mean_CASI_Co_L <- rowMeans(CPM[,CASI_Co_L], na.rm = T)

head(CPM[,CG_L])
head(CPM[,EUR_Cg_L])
head(CPM_S_Leaf)

# Average expression level in parents of CR & CO.
CPM_S_Leaf$Ave_CRCO <- rowMeans(CPM_S_Leaf[,c("mean_CR_L","mean_CO_L")])
CPM_S_Leaf$Ave_CGCO <- rowMeans(CPM_S_Leaf[,c("mean_CG_L","mean_CO_L")])
CPM_S_Leaf$Ave_CGCR <- rowMeans(CPM_S_Leaf[,c("mean_CG_L","mean_CR_L")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
CPM_S_Leaf$CPM_Scrco_Cbp_Cg <- (CPM_S_Leaf$mean_Cbp_Cg_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO # Subgenome Cg 

CPM_S_Leaf$CPM_Scrco_ASI_Cg <- (CPM_S_Leaf$mean_ASI_Cg_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_EUR_Cg <- (CPM_S_Leaf$mean_EUR_Cg_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_ME_Cg <- (CPM_S_Leaf$mean_ME_Cg_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_CASI_Cg <- (CPM_S_Leaf$mean_CASI_Cg_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO

CPM_S_Leaf$CPM_Scrco_Cbp_Co <- (CPM_S_Leaf$mean_Cbp_Co_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO # Subgenome Co

CPM_S_Leaf$CPM_Scrco_ASI_Co <- (CPM_S_Leaf$mean_ASI_Co_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_EUR_Co <- (CPM_S_Leaf$mean_EUR_Co_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_ME_Co <- (CPM_S_Leaf$mean_ME_Co_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO
CPM_S_Leaf$CPM_Scrco_CASI_Co <- (CPM_S_Leaf$mean_CASI_Co_L - CPM_S_Leaf$Ave_CRCO) / CPM_S_Leaf$Ave_CRCO

head(CPM_S_Leaf)
# in CG vs CO
CPM_S_Leaf$CPM_Scgco_Cbp_Cg <- (CPM_S_Leaf$mean_Cbp_Cg_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO  # Subgenome Cg 

CPM_S_Leaf$CPM_Scgco_ASI_Cg <- (CPM_S_Leaf$mean_ASI_Cg_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_EUR_Cg <- (CPM_S_Leaf$mean_EUR_Cg_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_ME_Cg <- (CPM_S_Leaf$mean_ME_Cg_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_CASI_Cg <- (CPM_S_Leaf$mean_CASI_Cg_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO

CPM_S_Leaf$CPM_Scgco_Cbp_Co <- (CPM_S_Leaf$mean_Cbp_Co_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO  # Subgenome Co

CPM_S_Leaf$CPM_Scgco_ASI_Co <- (CPM_S_Leaf$mean_ASI_Co_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_EUR_Co <- (CPM_S_Leaf$mean_EUR_Co_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_ME_Co <- (CPM_S_Leaf$mean_ME_Co_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO
CPM_S_Leaf$CPM_Scgco_CASI_Co <- (CPM_S_Leaf$mean_CASI_Co_L - CPM_S_Leaf$Ave_CGCO) / CPM_S_Leaf$Ave_CGCO

head(CPM_S_Leaf)
# in CG vs CR
CPM_S_Leaf$CPM_Scgcr_Cbp_Cg <- (CPM_S_Leaf$mean_Cbp_Cg_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR  # Subgenome Cg 

CPM_S_Leaf$CPM_Scgcr_ASI_Cg <- (CPM_S_Leaf$mean_ASI_Cg_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_EUR_Cg <- (CPM_S_Leaf$mean_EUR_Cg_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_ME_Cg <- (CPM_S_Leaf$mean_ME_Cg_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_CASI_Cg <- (CPM_S_Leaf$mean_CASI_Cg_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR

CPM_S_Leaf$CPM_Scgcr_Cbp_Co <- (CPM_S_Leaf$mean_Cbp_Co_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR  # Subgenome Co

CPM_S_Leaf$CPM_Scgcr_ASI_Co <- (CPM_S_Leaf$mean_ASI_Co_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_EUR_Co <- (CPM_S_Leaf$mean_EUR_Co_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_ME_Co <- (CPM_S_Leaf$mean_ME_Co_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR
CPM_S_Leaf$CPM_Scgcr_CASI_Co <- (CPM_S_Leaf$mean_CASI_Co_L - CPM_S_Leaf$Ave_CGCR) / CPM_S_Leaf$Ave_CGCR

head(CPM_S_Leaf)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_Cbp_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_Cbp_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_EUR_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_EUR_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ME_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ME_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_CASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_CASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_Cbp_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_Cbp_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ASI_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_EUR_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_EUR_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ME_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_ME_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_CASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CR_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scrco_CASI_Co * (-1)

# CG vs CO, CG > CO.

    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_Cbp_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_Cbp_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_EUR_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_EUR_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ME_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ME_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_CASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_CASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_Cbp_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_Cbp_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ASI_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_EUR_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_EUR_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ME_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_ME_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_CASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CO_L),]$CPM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.

    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_Cbp_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_Cbp_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_EUR_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_EUR_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ME_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ME_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_CASI_Cg <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_CASI_Cg * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_Cbp_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_Cbp_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ASI_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_EUR_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_EUR_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ME_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_ME_Co * (-1)
    CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_CASI_Co <- CPM_S_Leaf[(CPM_S_Leaf$mean_CG_L > CPM_S_Leaf$mean_CR_L),]$CPM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of CPM_S_Cbp_Cr and CPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
CPM_S_Leaf$delta_CPM_Scrco_Cbp <- abs(CPM_S_Leaf$CPM_Scrco_Cbp_Co) - abs(CPM_S_Leaf$CPM_Scrco_Cbp_Cg)
CPM_S_Leaf$delta_CPM_Scrco_ASI <- abs(CPM_S_Leaf$CPM_Scrco_ASI_Co) - abs(CPM_S_Leaf$CPM_Scrco_ASI_Cg)
CPM_S_Leaf$delta_CPM_Scrco_EUR <- abs(CPM_S_Leaf$CPM_Scrco_EUR_Co) - abs(CPM_S_Leaf$CPM_Scrco_EUR_Cg)
CPM_S_Leaf$delta_CPM_Scrco_ME <- abs(CPM_S_Leaf$CPM_Scrco_ME_Co) - abs(CPM_S_Leaf$CPM_Scrco_ME_Cg)
CPM_S_Leaf$delta_CPM_Scrco_CASI <- abs(CPM_S_Leaf$CPM_Scrco_CASI_Co) - abs(CPM_S_Leaf$CPM_Scrco_CASI_Cg)
head(CPM_S_Leaf)

# CG vs CO
CPM_S_Leaf$delta_CPM_Scgco_Cbp <- abs(CPM_S_Leaf$CPM_Scgco_Cbp_Co) - abs(CPM_S_Leaf$CPM_Scgco_Cbp_Cg)
CPM_S_Leaf$delta_CPM_Scgco_ASI <- abs(CPM_S_Leaf$CPM_Scgco_ASI_Co) - abs(CPM_S_Leaf$CPM_Scgco_ASI_Cg)
CPM_S_Leaf$delta_CPM_Scgco_EUR <- abs(CPM_S_Leaf$CPM_Scgco_EUR_Co) - abs(CPM_S_Leaf$CPM_Scgco_EUR_Cg)
CPM_S_Leaf$delta_CPM_Scgco_ME <- abs(CPM_S_Leaf$CPM_Scgco_ME_Co) - abs(CPM_S_Leaf$CPM_Scgco_ME_Cg)
CPM_S_Leaf$delta_CPM_Scgco_CASI <- abs(CPM_S_Leaf$CPM_Scgco_CASI_Co) - abs(CPM_S_Leaf$CPM_Scgco_CASI_Cg)
head(CPM_S_Leaf)

# CG vs CR
CPM_S_Leaf$delta_CPM_Scgcr_Cbp <- abs(CPM_S_Leaf$CPM_Scgcr_Cbp_Co) - abs(CPM_S_Leaf$CPM_Scgcr_Cbp_Cg)
CPM_S_Leaf$delta_CPM_Scgcr_ASI <- abs(CPM_S_Leaf$CPM_Scgcr_ASI_Co) - abs(CPM_S_Leaf$CPM_Scgcr_ASI_Cg)
CPM_S_Leaf$delta_CPM_Scgcr_EUR <- abs(CPM_S_Leaf$CPM_Scgcr_EUR_Co) - abs(CPM_S_Leaf$CPM_Scgcr_EUR_Cg)
CPM_S_Leaf$delta_CPM_Scgcr_ME <- abs(CPM_S_Leaf$CPM_Scgcr_ME_Co) - abs(CPM_S_Leaf$CPM_Scgcr_ME_Cg)
CPM_S_Leaf$delta_CPM_Scgcr_CASI <- abs(CPM_S_Leaf$CPM_Scgcr_CASI_Co) - abs(CPM_S_Leaf$CPM_Scgcr_CASI_Cg)
head(CPM_S_Leaf)
# load genes significantly expressed between Cr Co in flower.
DE.L.crco <- read.table("OutputData/DE_CrCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.crco)
dim(DE.L.crco)

DE.L.cgco <- read.table("OutputData/DE_CgCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgco)
dim(DE.L.cgco)

DE.L.cgcr <- read.table("OutputData/DE_CgCr_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgcr)
dim(DE.L.cgcr)

# keep genes with FDR ≤ 0.01
DE.L.crco.keep <- CPM_S_Leaf$Gene %in% row.names(DE.L.crco)
CPM_S_Leaf$Gene %in% row.names(DE.L.crco)
DE.L.cgco.keep <- CPM_S_Leaf$Gene %in% row.names(DE.L.cgco)
CPM_S_Leaf$Gene %in% row.names(DE.L.cgco)
DE.L.cgcr.keep <- CPM_S_Leaf$Gene %in% row.names(DE.L.cgcr)
CPM_S_Leaf$Gene %in% row.names(DE.L.cgcr)

CPM_Scrco_Leaf_FDR <- CPM_S_Leaf[DE.L.crco.keep,]
head(CPM_Scrco_Leaf_FDR)
dim(CPM_Scrco_Leaf_FDR)
CPM_Scgco_Leaf_FDR <- CPM_S_Leaf[DE.L.cgco.keep,]
head(CPM_Scgco_Leaf_FDR)
dim(CPM_Scgco_Leaf_FDR)
CPM_Scgcr_Leaf_FDR <- CPM_S_Leaf[DE.L.cgcr.keep,]
head(CPM_Scgcr_Leaf_FDR)
dim(CPM_Scgcr_Leaf_FDR)

#############################

# Create a new dataframe for Root
CPM_S_Root <- data.frame(Gene=rownames(CPM))
CPM_S_Root$mean_CG_R <- rowMeans(CPM[,CG_R], na.rm = T) # Average value of CG expression in Root
CPM_S_Root$mean_CR_R <- rowMeans(CPM[,CR_R], na.rm = T)
CPM_S_Root$mean_CO_R <- rowMeans(CPM[,CO_R], na.rm = T)

CPM_S_Root$mean_Cbp_Cg_R <- rowMeans(CPM[,Cbp_Cg_R], na.rm = T)
CPM_S_Root$mean_Cbp_Co_R <- rowMeans(CPM[,Cbp_Co_R], na.rm = T)

CPM_S_Root$mean_ASI_Cg_R <- rowMeans(CPM[,ASI_Cg_R], na.rm = T)
CPM_S_Root$mean_EUR_Cg_R <- rowMeans(CPM[,EUR_Cg_R], na.rm = T)
CPM_S_Root$mean_ME_Cg_R <- rowMeans(CPM[,ME_Cg_R], na.rm = T)
CPM_S_Root$mean_CASI_Cg_R <- rowMeans(CPM[,CASI_Cg_R], na.rm = T)

CPM_S_Root$mean_ASI_Co_R <- rowMeans(CPM[,ASI_Co_R], na.rm = T)
CPM_S_Root$mean_EUR_Co_R <- rowMeans(CPM[,EUR_Co_R], na.rm = T)
CPM_S_Root$mean_ME_Co_R <- rowMeans(CPM[,ME_Co_R], na.rm = T)
CPM_S_Root$mean_CASI_Co_R <- rowMeans(CPM[,CASI_Co_R], na.rm = T)

head(CPM[,CG_R])
head(CPM[,EUR_Cg_R])
head(CPM_S_Root)

# Average expression level in parents of CR & CO.
CPM_S_Root$Ave_CRCO <- rowMeans(CPM_S_Root[,c("mean_CR_R","mean_CO_R")])
CPM_S_Root$Ave_CGCO <- rowMeans(CPM_S_Root[,c("mean_CG_R","mean_CO_R")])
CPM_S_Root$Ave_CGCR <- rowMeans(CPM_S_Root[,c("mean_CG_R","mean_CR_R")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
CPM_S_Root$CPM_Scrco_Cbp_Cg <- (CPM_S_Root$mean_Cbp_Cg_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO # Subgenome Cg 

CPM_S_Root$CPM_Scrco_ASI_Cg <- (CPM_S_Root$mean_ASI_Cg_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_EUR_Cg <- (CPM_S_Root$mean_EUR_Cg_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_ME_Cg <- (CPM_S_Root$mean_ME_Cg_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_CASI_Cg <- (CPM_S_Root$mean_CASI_Cg_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO

CPM_S_Root$CPM_Scrco_Cbp_Co <- (CPM_S_Root$mean_Cbp_Co_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO # Subgenome Co

CPM_S_Root$CPM_Scrco_ASI_Co <- (CPM_S_Root$mean_ASI_Co_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_EUR_Co <- (CPM_S_Root$mean_EUR_Co_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_ME_Co <- (CPM_S_Root$mean_ME_Co_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO
CPM_S_Root$CPM_Scrco_CASI_Co <- (CPM_S_Root$mean_CASI_Co_R - CPM_S_Root$Ave_CRCO) / CPM_S_Root$Ave_CRCO

head(CPM_S_Root)
# in CG vs CO
CPM_S_Root$CPM_Scgco_Cbp_Cg <- (CPM_S_Root$mean_Cbp_Cg_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO  # Subgenome Cg 

CPM_S_Root$CPM_Scgco_ASI_Cg <- (CPM_S_Root$mean_ASI_Cg_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_EUR_Cg <- (CPM_S_Root$mean_EUR_Cg_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_ME_Cg <- (CPM_S_Root$mean_ME_Cg_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_CASI_Cg <- (CPM_S_Root$mean_CASI_Cg_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO

CPM_S_Root$CPM_Scgco_Cbp_Co <- (CPM_S_Root$mean_Cbp_Co_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO  # Subgenome Co

CPM_S_Root$CPM_Scgco_ASI_Co <- (CPM_S_Root$mean_ASI_Co_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_EUR_Co <- (CPM_S_Root$mean_EUR_Co_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_ME_Co <- (CPM_S_Root$mean_ME_Co_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO
CPM_S_Root$CPM_Scgco_CASI_Co <- (CPM_S_Root$mean_CASI_Co_R - CPM_S_Root$Ave_CGCO) / CPM_S_Root$Ave_CGCO

head(CPM_S_Root)
# in CG vs CR
CPM_S_Root$CPM_Scgcr_Cbp_Cg <- (CPM_S_Root$mean_Cbp_Cg_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR  # Subgenome Cg 

CPM_S_Root$CPM_Scgcr_ASI_Cg <- (CPM_S_Root$mean_ASI_Cg_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_EUR_Cg <- (CPM_S_Root$mean_EUR_Cg_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_ME_Cg <- (CPM_S_Root$mean_ME_Cg_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_CASI_Cg <- (CPM_S_Root$mean_CASI_Cg_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR

CPM_S_Root$CPM_Scgcr_Cbp_Co <- (CPM_S_Root$mean_Cbp_Co_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR  # Subgenome Co

CPM_S_Root$CPM_Scgcr_ASI_Co <- (CPM_S_Root$mean_ASI_Co_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_EUR_Co <- (CPM_S_Root$mean_EUR_Co_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_ME_Co <- (CPM_S_Root$mean_ME_Co_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR
CPM_S_Root$CPM_Scgcr_CASI_Co <- (CPM_S_Root$mean_CASI_Co_R - CPM_S_Root$Ave_CGCR) / CPM_S_Root$Ave_CGCR

head(CPM_S_Root)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_Cbp_Cg <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_Cbp_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_EUR_Cg <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_EUR_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ME_Cg <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ME_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_CASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_CASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_Cbp_Co <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_Cbp_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ASI_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_EUR_Co <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_EUR_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ME_Co <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_ME_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_CASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CR_R > CPM_S_Root$mean_CO_R),]$CPM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_Cbp_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_Cbp_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_EUR_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_EUR_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ME_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ME_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_CASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_CASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_Cbp_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_Cbp_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ASI_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_EUR_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_EUR_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ME_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_ME_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_CASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CO_R),]$CPM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_Cbp_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_Cbp_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_EUR_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_EUR_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ME_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ME_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_CASI_Cg <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_CASI_Cg * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_Cbp_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_Cbp_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ASI_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_EUR_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_EUR_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ME_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_ME_Co * (-1)
    CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_CASI_Co <- CPM_S_Root[(CPM_S_Root$mean_CG_R > CPM_S_Root$mean_CR_R),]$CPM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of CPM_S_Cbp_Cr and CPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
CPM_S_Root$delta_CPM_Scrco_Cbp <- abs(CPM_S_Root$CPM_Scrco_Cbp_Co) - abs(CPM_S_Root$CPM_Scrco_Cbp_Cg)
CPM_S_Root$delta_CPM_Scrco_ASI <- abs(CPM_S_Root$CPM_Scrco_ASI_Co) - abs(CPM_S_Root$CPM_Scrco_ASI_Cg)
CPM_S_Root$delta_CPM_Scrco_EUR <- abs(CPM_S_Root$CPM_Scrco_EUR_Co) - abs(CPM_S_Root$CPM_Scrco_EUR_Cg)
CPM_S_Root$delta_CPM_Scrco_ME <- abs(CPM_S_Root$CPM_Scrco_ME_Co) - abs(CPM_S_Root$CPM_Scrco_ME_Cg)
CPM_S_Root$delta_CPM_Scrco_CASI <- abs(CPM_S_Root$CPM_Scrco_CASI_Co) - abs(CPM_S_Root$CPM_Scrco_CASI_Cg)
head(CPM_S_Root)

# CG vs CO
CPM_S_Root$delta_CPM_Scgco_Cbp <- abs(CPM_S_Root$CPM_Scgco_Cbp_Co) - abs(CPM_S_Root$CPM_Scgco_Cbp_Cg)
CPM_S_Root$delta_CPM_Scgco_ASI <- abs(CPM_S_Root$CPM_Scgco_ASI_Co) - abs(CPM_S_Root$CPM_Scgco_ASI_Cg)
CPM_S_Root$delta_CPM_Scgco_EUR <- abs(CPM_S_Root$CPM_Scgco_EUR_Co) - abs(CPM_S_Root$CPM_Scgco_EUR_Cg)
CPM_S_Root$delta_CPM_Scgco_ME <- abs(CPM_S_Root$CPM_Scgco_ME_Co) - abs(CPM_S_Root$CPM_Scgco_ME_Cg)
CPM_S_Root$delta_CPM_Scgco_CASI <- abs(CPM_S_Root$CPM_Scgco_CASI_Co) - abs(CPM_S_Root$CPM_Scgco_CASI_Cg)
head(CPM_S_Root)

# CG vs CR
CPM_S_Root$delta_CPM_Scgcr_Cbp <- abs(CPM_S_Root$CPM_Scgcr_Cbp_Co) - abs(CPM_S_Root$CPM_Scgcr_Cbp_Cg)
CPM_S_Root$delta_CPM_Scgcr_ASI <- abs(CPM_S_Root$CPM_Scgcr_ASI_Co) - abs(CPM_S_Root$CPM_Scgcr_ASI_Cg)
CPM_S_Root$delta_CPM_Scgcr_EUR <- abs(CPM_S_Root$CPM_Scgcr_EUR_Co) - abs(CPM_S_Root$CPM_Scgcr_EUR_Cg)
CPM_S_Root$delta_CPM_Scgcr_ME <- abs(CPM_S_Root$CPM_Scgcr_ME_Co) - abs(CPM_S_Root$CPM_Scgcr_ME_Cg)
CPM_S_Root$delta_CPM_Scgcr_CASI <- abs(CPM_S_Root$CPM_Scgcr_CASI_Co) - abs(CPM_S_Root$CPM_Scgcr_CASI_Cg)
head(CPM_S_Root)
# load genes significantly expressed between Cr Co in flower.
DE.R.crco <- read.table("OutputData/DE_CrCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.crco)
dim(DE.R.crco)

DE.R.cgco <- read.table("OutputData/DE_CgCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgco)
dim(DE.R.cgco)

DE.R.cgcr <- read.table("OutputData/DE_CgCr_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgcr)
dim(DE.R.cgcr)

# keep genes with FDR ≤ 0.01
DE.R.crco.keep <- CPM_S_Root$Gene %in% row.names(DE.R.crco)
CPM_S_Root$Gene %in% row.names(DE.R.crco)
DE.R.cgco.keep <- CPM_S_Root$Gene %in% row.names(DE.R.cgco)
CPM_S_Root$Gene %in% row.names(DE.R.cgco)
DE.R.cgcr.keep <- CPM_S_Root$Gene %in% row.names(DE.R.cgcr)
CPM_S_Root$Gene %in% row.names(DE.R.cgcr)

CPM_Scrco_Root_FDR <- CPM_S_Root[DE.R.crco.keep,]
head(CPM_Scrco_Root_FDR)
dim(CPM_Scrco_Root_FDR)
CPM_Scgco_Root_FDR <- CPM_S_Root[DE.R.cgco.keep,]
head(CPM_Scgco_Root_FDR)
dim(CPM_Scgco_Root_FDR)
CPM_Scgcr_Root_FDR <- CPM_S_Root[DE.R.cgcr.keep,]
head(CPM_Scgcr_Root_FDR)
dim(CPM_Scgcr_Root_FDR)
########################################################################################################
#                                              logCPM                                                  #
########################################################################################################
# Similarity index
# logCPM
# Create a new dataframe for Flower
logCPM_S_Flower <- data.frame(Gene=rownames(logCPM))
logCPM_S_Flower$mean_CG_F <- rowMeans(logCPM[,CG_F], na.rm = T) # Average value of CG expression in Flower
logCPM_S_Flower$mean_CR_F <- rowMeans(logCPM[,CR_F], na.rm = T)
logCPM_S_Flower$mean_CO_F <- rowMeans(logCPM[,CO_F], na.rm = T)

logCPM_S_Flower$mean_Cbp_Cg_F <- rowMeans(logCPM[,Cbp_Cg_F], na.rm = T)
logCPM_S_Flower$mean_Cbp_Co_F <- rowMeans(logCPM[,Cbp_Co_F], na.rm = T)

logCPM_S_Flower$mean_ASI_Cg_F <- rowMeans(logCPM[,ASI_Cg_F], na.rm = T)
logCPM_S_Flower$mean_EUR_Cg_F <- rowMeans(logCPM[,EUR_Cg_F], na.rm = T)
logCPM_S_Flower$mean_ME_Cg_F <- rowMeans(logCPM[,ME_Cg_F], na.rm = T)
logCPM_S_Flower$mean_CASI_Cg_F <- rowMeans(logCPM[,CASI_Cg_F], na.rm = T)

logCPM_S_Flower$mean_ASI_Co_F <- rowMeans(logCPM[,ASI_Co_F], na.rm = T)
logCPM_S_Flower$mean_EUR_Co_F <- rowMeans(logCPM[,EUR_Co_F], na.rm = T)
logCPM_S_Flower$mean_ME_Co_F <- rowMeans(logCPM[,ME_Co_F], na.rm = T)
logCPM_S_Flower$mean_CASI_Co_F <- rowMeans(logCPM[,CASI_Co_F], na.rm = T)

head(logCPM[,CG_F])
head(logCPM[,EUR_Cg_F])
head(logCPM_S_Flower)

# Average expression level in parents of CR & CO.
logCPM_S_Flower$Ave_CRCO <- rowMeans(logCPM_S_Flower[,c("mean_CR_F","mean_CO_F")])
logCPM_S_Flower$Ave_CGCO <- rowMeans(logCPM_S_Flower[,c("mean_CG_F","mean_CO_F")])
logCPM_S_Flower$Ave_CGCR <- rowMeans(logCPM_S_Flower[,c("mean_CG_F","mean_CR_F")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
logCPM_S_Flower$logCPM_Scrco_Cbp_Cg <- (logCPM_S_Flower$mean_Cbp_Cg_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO # Subgenome Cg 

logCPM_S_Flower$logCPM_Scrco_ASI_Cg <- (logCPM_S_Flower$mean_ASI_Cg_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_EUR_Cg <- (logCPM_S_Flower$mean_EUR_Cg_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_ME_Cg <- (logCPM_S_Flower$mean_ME_Cg_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_CASI_Cg <- (logCPM_S_Flower$mean_CASI_Cg_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO

logCPM_S_Flower$logCPM_Scrco_Cbp_Co <- (logCPM_S_Flower$mean_Cbp_Co_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO # Subgenome Co

logCPM_S_Flower$logCPM_Scrco_ASI_Co <- (logCPM_S_Flower$mean_ASI_Co_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_EUR_Co <- (logCPM_S_Flower$mean_EUR_Co_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_ME_Co <- (logCPM_S_Flower$mean_ME_Co_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO
logCPM_S_Flower$logCPM_Scrco_CASI_Co <- (logCPM_S_Flower$mean_CASI_Co_F - logCPM_S_Flower$Ave_CRCO) / logCPM_S_Flower$Ave_CRCO

head(logCPM_S_Flower)
# in CG vs CO
logCPM_S_Flower$logCPM_Scgco_Cbp_Cg <- (logCPM_S_Flower$mean_Cbp_Cg_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO  # Subgenome Cg 

logCPM_S_Flower$logCPM_Scgco_ASI_Cg <- (logCPM_S_Flower$mean_ASI_Cg_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_EUR_Cg <- (logCPM_S_Flower$mean_EUR_Cg_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_ME_Cg <- (logCPM_S_Flower$mean_ME_Cg_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_CASI_Cg <- (logCPM_S_Flower$mean_CASI_Cg_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO

logCPM_S_Flower$logCPM_Scgco_Cbp_Co <- (logCPM_S_Flower$mean_Cbp_Co_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO  # Subgenome Co

logCPM_S_Flower$logCPM_Scgco_ASI_Co <- (logCPM_S_Flower$mean_ASI_Co_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_EUR_Co <- (logCPM_S_Flower$mean_EUR_Co_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_ME_Co <- (logCPM_S_Flower$mean_ME_Co_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO
logCPM_S_Flower$logCPM_Scgco_CASI_Co <- (logCPM_S_Flower$mean_CASI_Co_F - logCPM_S_Flower$Ave_CGCO) / logCPM_S_Flower$Ave_CGCO

head(logCPM_S_Flower)
# in CG vs CR
logCPM_S_Flower$logCPM_Scgcr_Cbp_Cg <- (logCPM_S_Flower$mean_Cbp_Cg_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR  # Subgenome Cg 

logCPM_S_Flower$logCPM_Scgcr_ASI_Cg <- (logCPM_S_Flower$mean_ASI_Cg_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_EUR_Cg <- (logCPM_S_Flower$mean_EUR_Cg_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_ME_Cg <- (logCPM_S_Flower$mean_ME_Cg_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_CASI_Cg <- (logCPM_S_Flower$mean_CASI_Cg_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR

logCPM_S_Flower$logCPM_Scgcr_Cbp_Co <- (logCPM_S_Flower$mean_Cbp_Co_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR  # Subgenome Co

logCPM_S_Flower$logCPM_Scgcr_ASI_Co <- (logCPM_S_Flower$mean_ASI_Co_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_EUR_Co <- (logCPM_S_Flower$mean_EUR_Co_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_ME_Co <- (logCPM_S_Flower$mean_ME_Co_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR
logCPM_S_Flower$logCPM_Scgcr_CASI_Co <- (logCPM_S_Flower$mean_CASI_Co_F - logCPM_S_Flower$Ave_CGCR) / logCPM_S_Flower$Ave_CGCR

head(logCPM_S_Flower)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_Cbp_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_Cbp_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_EUR_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_EUR_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ME_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ME_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_CASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_CASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_Cbp_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_Cbp_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ASI_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_EUR_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_EUR_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ME_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_ME_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_CASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CR_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_Cbp_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_Cbp_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_EUR_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_EUR_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ME_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ME_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_CASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_CASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_Cbp_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_Cbp_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ASI_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_EUR_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_EUR_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ME_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_ME_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_CASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CO_F),]$logCPM_Scgco_CASI_Co * (-1)



# CG vs CR, CG > CR.

logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_Cbp_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_Cbp_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_EUR_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_EUR_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ME_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ME_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_CASI_Cg <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_CASI_Cg * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_Cbp_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_Cbp_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ASI_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_EUR_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_EUR_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ME_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_ME_Co * (-1)
logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_CASI_Co <- logCPM_S_Flower[(logCPM_S_Flower$mean_CG_F > logCPM_S_Flower$mean_CR_F),]$logCPM_Scgcr_CASI_Co * (-1)

# Calculate delta S by comparing the absolute values of logCPM_S_Cbp_Cr and logCPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
logCPM_S_Flower$delta_logCPM_Scrco_Cbp <- abs(logCPM_S_Flower$logCPM_Scrco_Cbp_Cg) - abs(logCPM_S_Flower$logCPM_Scrco_Cbp_Co)
logCPM_S_Flower$delta_logCPM_Scrco_ASI <- abs(logCPM_S_Flower$logCPM_Scrco_ASI_Cg) - abs(logCPM_S_Flower$logCPM_Scrco_ASI_Co)
logCPM_S_Flower$delta_logCPM_Scrco_EUR <- abs(logCPM_S_Flower$logCPM_Scrco_EUR_Cg) - abs(logCPM_S_Flower$logCPM_Scrco_EUR_Co)
logCPM_S_Flower$delta_logCPM_Scrco_ME <- abs(logCPM_S_Flower$logCPM_Scrco_ME_Cg) - abs(logCPM_S_Flower$logCPM_Scrco_ME_Co)
logCPM_S_Flower$delta_logCPM_Scrco_CASI <- abs(logCPM_S_Flower$logCPM_Scrco_CASI_Cg) - abs(logCPM_S_Flower$logCPM_Scrco_CASI_Co)
head(logCPM_S_Flower)

# CG vs CO
logCPM_S_Flower$delta_logCPM_Scgco_Cbp <- abs(logCPM_S_Flower$logCPM_Scgco_Cbp_Cg) - abs(logCPM_S_Flower$logCPM_Scgco_Cbp_Co)
logCPM_S_Flower$delta_logCPM_Scgco_ASI <- abs(logCPM_S_Flower$logCPM_Scgco_ASI_Cg) - abs(logCPM_S_Flower$logCPM_Scgco_ASI_Co)
logCPM_S_Flower$delta_logCPM_Scgco_EUR <- abs(logCPM_S_Flower$logCPM_Scgco_EUR_Cg) - abs(logCPM_S_Flower$logCPM_Scgco_EUR_Co)
logCPM_S_Flower$delta_logCPM_Scgco_ME <- abs(logCPM_S_Flower$logCPM_Scgco_ME_Cg) - abs(logCPM_S_Flower$logCPM_Scgco_ME_Co)
logCPM_S_Flower$delta_logCPM_Scgco_CASI <- abs(logCPM_S_Flower$logCPM_Scgco_CASI_Cg) - abs(logCPM_S_Flower$logCPM_Scgco_CASI_Co)
head(logCPM_S_Flower)

# CG vs CR
logCPM_S_Flower$delta_logCPM_Scgcr_Cbp <- abs(logCPM_S_Flower$logCPM_Scgcr_Cbp_Cg) - abs(logCPM_S_Flower$logCPM_Scgcr_Cbp_Co)
logCPM_S_Flower$delta_logCPM_Scgcr_ASI <- abs(logCPM_S_Flower$logCPM_Scgcr_ASI_Cg) - abs(logCPM_S_Flower$logCPM_Scgcr_ASI_Co)
logCPM_S_Flower$delta_logCPM_Scgcr_EUR <- abs(logCPM_S_Flower$logCPM_Scgcr_EUR_Cg) - abs(logCPM_S_Flower$logCPM_Scgcr_EUR_Co)
logCPM_S_Flower$delta_logCPM_Scgcr_ME <- abs(logCPM_S_Flower$logCPM_Scgcr_ME_Cg) - abs(logCPM_S_Flower$logCPM_Scgcr_ME_Co)
logCPM_S_Flower$delta_logCPM_Scgcr_CASI <- abs(logCPM_S_Flower$logCPM_Scgcr_CASI_Cg) - abs(logCPM_S_Flower$logCPM_Scgcr_CASI_Co)
head(logCPM_S_Flower)
# load genes significantly expressed between Cr Co in flower.
DE.F.crco <- read.table("OutputData/DE_CrCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.crco)
dim(DE.F.crco)

DE.F.cgco <- read.table("OutputData/DE_CgCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgco)
dim(DE.F.cgco)

DE.F.cgcr <- read.table("OutputData/DE_CgCr_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgcr)
dim(DE.F.cgcr)

# keep genes with FDR ≤ 0.01
DE.F.crco.keep <- logCPM_S_Flower$Gene %in% row.names(DE.F.crco)
logCPM_S_Flower$Gene %in% row.names(DE.F.crco)
DE.F.cgco.keep <- logCPM_S_Flower$Gene %in% row.names(DE.F.cgco)
logCPM_S_Flower$Gene %in% row.names(DE.F.cgco)
DE.F.cgcr.keep <- logCPM_S_Flower$Gene %in% row.names(DE.F.cgcr)
logCPM_S_Flower$Gene %in% row.names(DE.F.cgcr)

logCPM_Scrco_Flower_FDR <- logCPM_S_Flower[DE.F.crco.keep,]
head(logCPM_Scrco_Flower_FDR)
dim(logCPM_Scrco_Flower_FDR)
logCPM_Scgco_Flower_FDR <- logCPM_S_Flower[DE.F.cgco.keep,]
head(logCPM_Scgco_Flower_FDR)
dim(logCPM_Scgco_Flower_FDR)
logCPM_Scgcr_Flower_FDR <- logCPM_S_Flower[DE.F.cgcr.keep,]
head(logCPM_Scgcr_Flower_FDR)
dim(logCPM_Scgcr_Flower_FDR)

##################################

# Create a new dataframe for Leaf
logCPM_S_Leaf <- data.frame(Gene=rownames(logCPM))
logCPM_S_Leaf$mean_CG_L <- rowMeans(logCPM[,CG_L], na.rm = T) # Average value of CG expression in Leaf
logCPM_S_Leaf$mean_CR_L <- rowMeans(logCPM[,CR_L], na.rm = T)
logCPM_S_Leaf$mean_CO_L <- rowMeans(logCPM[,CO_L], na.rm = T)

logCPM_S_Leaf$mean_Cbp_Cg_L <- rowMeans(logCPM[,Cbp_Cg_L], na.rm = T)
logCPM_S_Leaf$mean_Cbp_Co_L <- rowMeans(logCPM[,Cbp_Co_L], na.rm = T)

logCPM_S_Leaf$mean_ASI_Cg_L <- rowMeans(logCPM[,ASI_Cg_L], na.rm = T)
logCPM_S_Leaf$mean_EUR_Cg_L <- rowMeans(logCPM[,EUR_Cg_L], na.rm = T)
logCPM_S_Leaf$mean_ME_Cg_L <- rowMeans(logCPM[,ME_Cg_L], na.rm = T)
logCPM_S_Leaf$mean_CASI_Cg_L <- rowMeans(logCPM[,CASI_Cg_L], na.rm = T)

logCPM_S_Leaf$mean_ASI_Co_L <- rowMeans(logCPM[,ASI_Co_L], na.rm = T)
logCPM_S_Leaf$mean_EUR_Co_L <- rowMeans(logCPM[,EUR_Co_L], na.rm = T)
logCPM_S_Leaf$mean_ME_Co_L <- rowMeans(logCPM[,ME_Co_L], na.rm = T)
logCPM_S_Leaf$mean_CASI_Co_L <- rowMeans(logCPM[,CASI_Co_L], na.rm = T)

head(logCPM[,CG_L])
head(logCPM[,EUR_Cg_L])
head(logCPM_S_Leaf)

# Average expression level in parents of CR & CO.
logCPM_S_Leaf$Ave_CRCO <- rowMeans(logCPM_S_Leaf[,c("mean_CR_L","mean_CO_L")])
logCPM_S_Leaf$Ave_CGCO <- rowMeans(logCPM_S_Leaf[,c("mean_CG_L","mean_CO_L")])
logCPM_S_Leaf$Ave_CGCR <- rowMeans(logCPM_S_Leaf[,c("mean_CG_L","mean_CR_L")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
logCPM_S_Leaf$logCPM_Scrco_Cbp_Cg <- (logCPM_S_Leaf$mean_Cbp_Cg_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO # Subgenome Cg 

logCPM_S_Leaf$logCPM_Scrco_ASI_Cg <- (logCPM_S_Leaf$mean_ASI_Cg_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_EUR_Cg <- (logCPM_S_Leaf$mean_EUR_Cg_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_ME_Cg <- (logCPM_S_Leaf$mean_ME_Cg_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_CASI_Cg <- (logCPM_S_Leaf$mean_CASI_Cg_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO

logCPM_S_Leaf$logCPM_Scrco_Cbp_Co <- (logCPM_S_Leaf$mean_Cbp_Co_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO # Subgenome Co

logCPM_S_Leaf$logCPM_Scrco_ASI_Co <- (logCPM_S_Leaf$mean_ASI_Co_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_EUR_Co <- (logCPM_S_Leaf$mean_EUR_Co_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_ME_Co <- (logCPM_S_Leaf$mean_ME_Co_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO
logCPM_S_Leaf$logCPM_Scrco_CASI_Co <- (logCPM_S_Leaf$mean_CASI_Co_L - logCPM_S_Leaf$Ave_CRCO) / logCPM_S_Leaf$Ave_CRCO

head(logCPM_S_Leaf)
# in CG vs CO
logCPM_S_Leaf$logCPM_Scgco_Cbp_Cg <- (logCPM_S_Leaf$mean_Cbp_Cg_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO  # Subgenome Cg 

logCPM_S_Leaf$logCPM_Scgco_ASI_Cg <- (logCPM_S_Leaf$mean_ASI_Cg_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_EUR_Cg <- (logCPM_S_Leaf$mean_EUR_Cg_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_ME_Cg <- (logCPM_S_Leaf$mean_ME_Cg_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_CASI_Cg <- (logCPM_S_Leaf$mean_CASI_Cg_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO

logCPM_S_Leaf$logCPM_Scgco_Cbp_Co <- (logCPM_S_Leaf$mean_Cbp_Co_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO  # Subgenome Co

logCPM_S_Leaf$logCPM_Scgco_ASI_Co <- (logCPM_S_Leaf$mean_ASI_Co_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_EUR_Co <- (logCPM_S_Leaf$mean_EUR_Co_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_ME_Co <- (logCPM_S_Leaf$mean_ME_Co_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO
logCPM_S_Leaf$logCPM_Scgco_CASI_Co <- (logCPM_S_Leaf$mean_CASI_Co_L - logCPM_S_Leaf$Ave_CGCO) / logCPM_S_Leaf$Ave_CGCO

head(logCPM_S_Leaf)
# in CG vs CR
logCPM_S_Leaf$logCPM_Scgcr_Cbp_Cg <- (logCPM_S_Leaf$mean_Cbp_Cg_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR  # Subgenome Cg 

logCPM_S_Leaf$logCPM_Scgcr_ASI_Cg <- (logCPM_S_Leaf$mean_ASI_Cg_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_EUR_Cg <- (logCPM_S_Leaf$mean_EUR_Cg_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_ME_Cg <- (logCPM_S_Leaf$mean_ME_Cg_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_CASI_Cg <- (logCPM_S_Leaf$mean_CASI_Cg_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR

logCPM_S_Leaf$logCPM_Scgcr_Cbp_Co <- (logCPM_S_Leaf$mean_Cbp_Co_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR  # Subgenome Co

logCPM_S_Leaf$logCPM_Scgcr_ASI_Co <- (logCPM_S_Leaf$mean_ASI_Co_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_EUR_Co <- (logCPM_S_Leaf$mean_EUR_Co_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_ME_Co <- (logCPM_S_Leaf$mean_ME_Co_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR
logCPM_S_Leaf$logCPM_Scgcr_CASI_Co <- (logCPM_S_Leaf$mean_CASI_Co_L - logCPM_S_Leaf$Ave_CGCR) / logCPM_S_Leaf$Ave_CGCR

head(logCPM_S_Leaf)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_Cbp_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_Cbp_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_EUR_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_EUR_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ME_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ME_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_CASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_CASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_Cbp_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_Cbp_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ASI_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_EUR_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_EUR_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ME_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_ME_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_CASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CR_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scrco_CASI_Co * (-1)

# CG vs CO, CG > CO.

logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_Cbp_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_Cbp_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_EUR_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_EUR_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ME_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ME_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_CASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_CASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_Cbp_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_Cbp_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ASI_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_EUR_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_EUR_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ME_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_ME_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_CASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CO_L),]$logCPM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.

logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_Cbp_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_Cbp_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_EUR_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_EUR_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ME_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ME_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_CASI_Cg <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_CASI_Cg * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_Cbp_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_Cbp_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ASI_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_EUR_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_EUR_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ME_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_ME_Co * (-1)
logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_CASI_Co <- logCPM_S_Leaf[(logCPM_S_Leaf$mean_CG_L > logCPM_S_Leaf$mean_CR_L),]$logCPM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of logCPM_S_Cbp_Cr and logCPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
logCPM_S_Leaf$delta_logCPM_Scrco_Cbp <- abs(logCPM_S_Leaf$logCPM_Scrco_Cbp_Co) - abs(logCPM_S_Leaf$logCPM_Scrco_Cbp_Cg)
logCPM_S_Leaf$delta_logCPM_Scrco_ASI <- abs(logCPM_S_Leaf$logCPM_Scrco_ASI_Co) - abs(logCPM_S_Leaf$logCPM_Scrco_ASI_Cg)
logCPM_S_Leaf$delta_logCPM_Scrco_EUR <- abs(logCPM_S_Leaf$logCPM_Scrco_EUR_Co) - abs(logCPM_S_Leaf$logCPM_Scrco_EUR_Cg)
logCPM_S_Leaf$delta_logCPM_Scrco_ME <- abs(logCPM_S_Leaf$logCPM_Scrco_ME_Co) - abs(logCPM_S_Leaf$logCPM_Scrco_ME_Cg)
logCPM_S_Leaf$delta_logCPM_Scrco_CASI <- abs(logCPM_S_Leaf$logCPM_Scrco_CASI_Co) - abs(logCPM_S_Leaf$logCPM_Scrco_CASI_Cg)
head(logCPM_S_Leaf)

# CG vs CO
logCPM_S_Leaf$delta_logCPM_Scgco_Cbp <- abs(logCPM_S_Leaf$logCPM_Scgco_Cbp_Co) - abs(logCPM_S_Leaf$logCPM_Scgco_Cbp_Cg)
logCPM_S_Leaf$delta_logCPM_Scgco_ASI <- abs(logCPM_S_Leaf$logCPM_Scgco_ASI_Co) - abs(logCPM_S_Leaf$logCPM_Scgco_ASI_Cg)
logCPM_S_Leaf$delta_logCPM_Scgco_EUR <- abs(logCPM_S_Leaf$logCPM_Scgco_EUR_Co) - abs(logCPM_S_Leaf$logCPM_Scgco_EUR_Cg)
logCPM_S_Leaf$delta_logCPM_Scgco_ME <- abs(logCPM_S_Leaf$logCPM_Scgco_ME_Co) - abs(logCPM_S_Leaf$logCPM_Scgco_ME_Cg)
logCPM_S_Leaf$delta_logCPM_Scgco_CASI <- abs(logCPM_S_Leaf$logCPM_Scgco_CASI_Co) - abs(logCPM_S_Leaf$logCPM_Scgco_CASI_Cg)
head(logCPM_S_Leaf)

# CG vs CR
logCPM_S_Leaf$delta_logCPM_Scgcr_Cbp <- abs(logCPM_S_Leaf$logCPM_Scgcr_Cbp_Co) - abs(logCPM_S_Leaf$logCPM_Scgcr_Cbp_Cg)
logCPM_S_Leaf$delta_logCPM_Scgcr_ASI <- abs(logCPM_S_Leaf$logCPM_Scgcr_ASI_Co) - abs(logCPM_S_Leaf$logCPM_Scgcr_ASI_Cg)
logCPM_S_Leaf$delta_logCPM_Scgcr_EUR <- abs(logCPM_S_Leaf$logCPM_Scgcr_EUR_Co) - abs(logCPM_S_Leaf$logCPM_Scgcr_EUR_Cg)
logCPM_S_Leaf$delta_logCPM_Scgcr_ME <- abs(logCPM_S_Leaf$logCPM_Scgcr_ME_Co) - abs(logCPM_S_Leaf$logCPM_Scgcr_ME_Cg)
logCPM_S_Leaf$delta_logCPM_Scgcr_CASI <- abs(logCPM_S_Leaf$logCPM_Scgcr_CASI_Co) - abs(logCPM_S_Leaf$logCPM_Scgcr_CASI_Cg)
head(logCPM_S_Leaf)
# load genes significantly expressed between Cr Co in flower.
DE.L.crco <- read.table("OutputData/DE_CrCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.crco)
dim(DE.L.crco)

DE.L.cgco <- read.table("OutputData/DE_CgCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgco)
dim(DE.L.cgco)

DE.L.cgcr <- read.table("OutputData/DE_CgCr_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgcr)
dim(DE.L.cgcr)

# keep genes with FDR ≤ 0.01
DE.L.crco.keep <- logCPM_S_Leaf$Gene %in% row.names(DE.L.crco)
logCPM_S_Leaf$Gene %in% row.names(DE.L.crco)
DE.L.cgco.keep <- logCPM_S_Leaf$Gene %in% row.names(DE.L.cgco)
logCPM_S_Leaf$Gene %in% row.names(DE.L.cgco)
DE.L.cgcr.keep <- logCPM_S_Leaf$Gene %in% row.names(DE.L.cgcr)
logCPM_S_Leaf$Gene %in% row.names(DE.L.cgcr)

logCPM_Scrco_Leaf_FDR <- logCPM_S_Leaf[DE.L.crco.keep,]
head(logCPM_Scrco_Leaf_FDR)
dim(logCPM_Scrco_Leaf_FDR)
logCPM_Scgco_Leaf_FDR <- logCPM_S_Leaf[DE.L.cgco.keep,]
head(logCPM_Scgco_Leaf_FDR)
dim(logCPM_Scgco_Leaf_FDR)
logCPM_Scgcr_Leaf_FDR <- logCPM_S_Leaf[DE.L.cgcr.keep,]
head(logCPM_Scgcr_Leaf_FDR)
dim(logCPM_Scgcr_Leaf_FDR)

#############################

# Create a new dataframe for Root
logCPM_S_Root <- data.frame(Gene=rownames(logCPM))
logCPM_S_Root$mean_CG_R <- rowMeans(logCPM[,CG_R], na.rm = T) # Average value of CG expression in Root
logCPM_S_Root$mean_CR_R <- rowMeans(logCPM[,CR_R], na.rm = T)
logCPM_S_Root$mean_CO_R <- rowMeans(logCPM[,CO_R], na.rm = T)

logCPM_S_Root$mean_Cbp_Cg_R <- rowMeans(logCPM[,Cbp_Cg_R], na.rm = T)
logCPM_S_Root$mean_Cbp_Co_R <- rowMeans(logCPM[,Cbp_Co_R], na.rm = T)

logCPM_S_Root$mean_ASI_Cg_R <- rowMeans(logCPM[,ASI_Cg_R], na.rm = T)
logCPM_S_Root$mean_EUR_Cg_R <- rowMeans(logCPM[,EUR_Cg_R], na.rm = T)
logCPM_S_Root$mean_ME_Cg_R <- rowMeans(logCPM[,ME_Cg_R], na.rm = T)
logCPM_S_Root$mean_CASI_Cg_R <- rowMeans(logCPM[,CASI_Cg_R], na.rm = T)

logCPM_S_Root$mean_ASI_Co_R <- rowMeans(logCPM[,ASI_Co_R], na.rm = T)
logCPM_S_Root$mean_EUR_Co_R <- rowMeans(logCPM[,EUR_Co_R], na.rm = T)
logCPM_S_Root$mean_ME_Co_R <- rowMeans(logCPM[,ME_Co_R], na.rm = T)
logCPM_S_Root$mean_CASI_Co_R <- rowMeans(logCPM[,CASI_Co_R], na.rm = T)

head(logCPM[,CG_R])
head(logCPM[,EUR_Cg_R])
head(logCPM_S_Root)

# Average expression level in parents of CR & CO.
logCPM_S_Root$Ave_CRCO <- rowMeans(logCPM_S_Root[,c("mean_CR_R","mean_CO_R")])
logCPM_S_Root$Ave_CGCO <- rowMeans(logCPM_S_Root[,c("mean_CG_R","mean_CO_R")])
logCPM_S_Root$Ave_CGCR <- rowMeans(logCPM_S_Root[,c("mean_CG_R","mean_CR_R")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
logCPM_S_Root$logCPM_Scrco_Cbp_Cg <- (logCPM_S_Root$mean_Cbp_Cg_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO # Subgenome Cg 

logCPM_S_Root$logCPM_Scrco_ASI_Cg <- (logCPM_S_Root$mean_ASI_Cg_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_EUR_Cg <- (logCPM_S_Root$mean_EUR_Cg_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_ME_Cg <- (logCPM_S_Root$mean_ME_Cg_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_CASI_Cg <- (logCPM_S_Root$mean_CASI_Cg_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO

logCPM_S_Root$logCPM_Scrco_Cbp_Co <- (logCPM_S_Root$mean_Cbp_Co_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO # Subgenome Co

logCPM_S_Root$logCPM_Scrco_ASI_Co <- (logCPM_S_Root$mean_ASI_Co_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_EUR_Co <- (logCPM_S_Root$mean_EUR_Co_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_ME_Co <- (logCPM_S_Root$mean_ME_Co_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO
logCPM_S_Root$logCPM_Scrco_CASI_Co <- (logCPM_S_Root$mean_CASI_Co_R - logCPM_S_Root$Ave_CRCO) / logCPM_S_Root$Ave_CRCO

head(logCPM_S_Root)
# in CG vs CO
logCPM_S_Root$logCPM_Scgco_Cbp_Cg <- (logCPM_S_Root$mean_Cbp_Cg_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO  # Subgenome Cg 

logCPM_S_Root$logCPM_Scgco_ASI_Cg <- (logCPM_S_Root$mean_ASI_Cg_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_EUR_Cg <- (logCPM_S_Root$mean_EUR_Cg_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_ME_Cg <- (logCPM_S_Root$mean_ME_Cg_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_CASI_Cg <- (logCPM_S_Root$mean_CASI_Cg_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO

logCPM_S_Root$logCPM_Scgco_Cbp_Co <- (logCPM_S_Root$mean_Cbp_Co_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO  # Subgenome Co

logCPM_S_Root$logCPM_Scgco_ASI_Co <- (logCPM_S_Root$mean_ASI_Co_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_EUR_Co <- (logCPM_S_Root$mean_EUR_Co_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_ME_Co <- (logCPM_S_Root$mean_ME_Co_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO
logCPM_S_Root$logCPM_Scgco_CASI_Co <- (logCPM_S_Root$mean_CASI_Co_R - logCPM_S_Root$Ave_CGCO) / logCPM_S_Root$Ave_CGCO

head(logCPM_S_Root)
# in CG vs CR
logCPM_S_Root$logCPM_Scgcr_Cbp_Cg <- (logCPM_S_Root$mean_Cbp_Cg_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR  # Subgenome Cg 

logCPM_S_Root$logCPM_Scgcr_ASI_Cg <- (logCPM_S_Root$mean_ASI_Cg_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_EUR_Cg <- (logCPM_S_Root$mean_EUR_Cg_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_ME_Cg <- (logCPM_S_Root$mean_ME_Cg_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_CASI_Cg <- (logCPM_S_Root$mean_CASI_Cg_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR

logCPM_S_Root$logCPM_Scgcr_Cbp_Co <- (logCPM_S_Root$mean_Cbp_Co_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR  # Subgenome Co

logCPM_S_Root$logCPM_Scgcr_ASI_Co <- (logCPM_S_Root$mean_ASI_Co_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_EUR_Co <- (logCPM_S_Root$mean_EUR_Co_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_ME_Co <- (logCPM_S_Root$mean_ME_Co_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR
logCPM_S_Root$logCPM_Scgcr_CASI_Co <- (logCPM_S_Root$mean_CASI_Co_R - logCPM_S_Root$Ave_CGCR) / logCPM_S_Root$Ave_CGCR

head(logCPM_S_Root)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_Cbp_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_Cbp_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_EUR_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_EUR_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ME_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ME_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_CASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_CASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_Cbp_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_Cbp_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ASI_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_EUR_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_EUR_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ME_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_ME_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_CASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CR_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_Cbp_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_Cbp_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_EUR_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_EUR_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ME_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ME_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_CASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_CASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_Cbp_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_Cbp_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ASI_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_EUR_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_EUR_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ME_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_ME_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_CASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CO_R),]$logCPM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_Cbp_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_Cbp_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_EUR_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_EUR_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ME_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ME_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_CASI_Cg <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_CASI_Cg * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_Cbp_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_Cbp_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ASI_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_EUR_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_EUR_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ME_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_ME_Co * (-1)
logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_CASI_Co <- logCPM_S_Root[(logCPM_S_Root$mean_CG_R > logCPM_S_Root$mean_CR_R),]$logCPM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of logCPM_S_Cbp_Cr and logCPM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
logCPM_S_Root$delta_logCPM_Scrco_Cbp <- abs(logCPM_S_Root$logCPM_Scrco_Cbp_Co) - abs(logCPM_S_Root$logCPM_Scrco_Cbp_Cg)
logCPM_S_Root$delta_logCPM_Scrco_ASI <- abs(logCPM_S_Root$logCPM_Scrco_ASI_Co) - abs(logCPM_S_Root$logCPM_Scrco_ASI_Cg)
logCPM_S_Root$delta_logCPM_Scrco_EUR <- abs(logCPM_S_Root$logCPM_Scrco_EUR_Co) - abs(logCPM_S_Root$logCPM_Scrco_EUR_Cg)
logCPM_S_Root$delta_logCPM_Scrco_ME <- abs(logCPM_S_Root$logCPM_Scrco_ME_Co) - abs(logCPM_S_Root$logCPM_Scrco_ME_Cg)
logCPM_S_Root$delta_logCPM_Scrco_CASI <- abs(logCPM_S_Root$logCPM_Scrco_CASI_Co) - abs(logCPM_S_Root$logCPM_Scrco_CASI_Cg)
head(logCPM_S_Root)

# CG vs CO
logCPM_S_Root$delta_logCPM_Scgco_Cbp <- abs(logCPM_S_Root$logCPM_Scgco_Cbp_Co) - abs(logCPM_S_Root$logCPM_Scgco_Cbp_Cg)
logCPM_S_Root$delta_logCPM_Scgco_ASI <- abs(logCPM_S_Root$logCPM_Scgco_ASI_Co) - abs(logCPM_S_Root$logCPM_Scgco_ASI_Cg)
logCPM_S_Root$delta_logCPM_Scgco_EUR <- abs(logCPM_S_Root$logCPM_Scgco_EUR_Co) - abs(logCPM_S_Root$logCPM_Scgco_EUR_Cg)
logCPM_S_Root$delta_logCPM_Scgco_ME <- abs(logCPM_S_Root$logCPM_Scgco_ME_Co) - abs(logCPM_S_Root$logCPM_Scgco_ME_Cg)
logCPM_S_Root$delta_logCPM_Scgco_CASI <- abs(logCPM_S_Root$logCPM_Scgco_CASI_Co) - abs(logCPM_S_Root$logCPM_Scgco_CASI_Cg)
head(logCPM_S_Root)

# CG vs CR
logCPM_S_Root$delta_logCPM_Scgcr_Cbp <- abs(logCPM_S_Root$logCPM_Scgcr_Cbp_Co) - abs(logCPM_S_Root$logCPM_Scgcr_Cbp_Cg)
logCPM_S_Root$delta_logCPM_Scgcr_ASI <- abs(logCPM_S_Root$logCPM_Scgcr_ASI_Co) - abs(logCPM_S_Root$logCPM_Scgcr_ASI_Cg)
logCPM_S_Root$delta_logCPM_Scgcr_EUR <- abs(logCPM_S_Root$logCPM_Scgcr_EUR_Co) - abs(logCPM_S_Root$logCPM_Scgcr_EUR_Cg)
logCPM_S_Root$delta_logCPM_Scgcr_ME <- abs(logCPM_S_Root$logCPM_Scgcr_ME_Co) - abs(logCPM_S_Root$logCPM_Scgcr_ME_Cg)
logCPM_S_Root$delta_logCPM_Scgcr_CASI <- abs(logCPM_S_Root$logCPM_Scgcr_CASI_Co) - abs(logCPM_S_Root$logCPM_Scgcr_CASI_Cg)
head(logCPM_S_Root)
# load genes significantly expressed between Cr Co in flower.
DE.R.crco <- read.table("OutputData/DE_CrCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.crco)
dim(DE.R.crco)

DE.R.cgco <- read.table("OutputData/DE_CgCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgco)
dim(DE.R.cgco)

DE.R.cgcr <- read.table("OutputData/DE_CgCr_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgcr)
dim(DE.R.cgcr)

# keep genes with FDR ≤ 0.01
DE.R.crco.keep <- logCPM_S_Root$Gene %in% row.names(DE.R.crco)
logCPM_S_Root$Gene %in% row.names(DE.R.crco)
DE.R.cgco.keep <- logCPM_S_Root$Gene %in% row.names(DE.R.cgco)
logCPM_S_Root$Gene %in% row.names(DE.R.cgco)
DE.R.cgcr.keep <- logCPM_S_Root$Gene %in% row.names(DE.R.cgcr)
logCPM_S_Root$Gene %in% row.names(DE.R.cgcr)

logCPM_Scrco_Root_FDR <- logCPM_S_Root[DE.R.crco.keep,]
head(logCPM_Scrco_Root_FDR)
dim(logCPM_Scrco_Root_FDR)
logCPM_Scgco_Root_FDR <- logCPM_S_Root[DE.R.cgco.keep,]
head(logCPM_Scgco_Root_FDR)
dim(logCPM_Scgco_Root_FDR)
logCPM_Scgcr_Root_FDR <- logCPM_S_Root[DE.R.cgcr.keep,]
head(logCPM_Scgcr_Root_FDR)
dim(logCPM_Scgcr_Root_FDR)
########################################################################################################
#                                               TMM                                                    #
########################################################################################################
# Similarity index
# TMM
# Create a new dataframe for Flower
TMM_S_Flower <- data.frame(Gene=rownames(TMM))
TMM_S_Flower$mean_CG_F <- rowMeans(TMM[,CG_F], na.rm = T) # Average value of CG expression in Flower
TMM_S_Flower$mean_CR_F <- rowMeans(TMM[,CR_F], na.rm = T)
TMM_S_Flower$mean_CO_F <- rowMeans(TMM[,CO_F], na.rm = T)

TMM_S_Flower$mean_Cbp_Cg_F <- rowMeans(TMM[,Cbp_Cg_F], na.rm = T)
TMM_S_Flower$mean_Cbp_Co_F <- rowMeans(TMM[,Cbp_Co_F], na.rm = T)

TMM_S_Flower$mean_ASI_Cg_F <- rowMeans(TMM[,ASI_Cg_F], na.rm = T)
TMM_S_Flower$mean_EUR_Cg_F <- rowMeans(TMM[,EUR_Cg_F], na.rm = T)
TMM_S_Flower$mean_ME_Cg_F <- rowMeans(TMM[,ME_Cg_F], na.rm = T)
TMM_S_Flower$mean_CASI_Cg_F <- rowMeans(TMM[,CASI_Cg_F], na.rm = T)

TMM_S_Flower$mean_ASI_Co_F <- rowMeans(TMM[,ASI_Co_F], na.rm = T)
TMM_S_Flower$mean_EUR_Co_F <- rowMeans(TMM[,EUR_Co_F], na.rm = T)
TMM_S_Flower$mean_ME_Co_F <- rowMeans(TMM[,ME_Co_F], na.rm = T)
TMM_S_Flower$mean_CASI_Co_F <- rowMeans(TMM[,CASI_Co_F], na.rm = T)

head(TMM[,CG_F])
head(TMM[,EUR_Cg_F])
head(TMM_S_Flower)

# Average expression level in parents of CR & CO.
TMM_S_Flower$Ave_CRCO <- rowMeans(TMM_S_Flower[,c("mean_CR_F","mean_CO_F")])
TMM_S_Flower$Ave_CGCO <- rowMeans(TMM_S_Flower[,c("mean_CG_F","mean_CO_F")])
TMM_S_Flower$Ave_CGCR <- rowMeans(TMM_S_Flower[,c("mean_CG_F","mean_CR_F")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
TMM_S_Flower$TMM_Scrco_Cbp_Cg <- (TMM_S_Flower$mean_Cbp_Cg_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO # Subgenome Cg 

TMM_S_Flower$TMM_Scrco_ASI_Cg <- (TMM_S_Flower$mean_ASI_Cg_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_EUR_Cg <- (TMM_S_Flower$mean_EUR_Cg_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_ME_Cg <- (TMM_S_Flower$mean_ME_Cg_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_CASI_Cg <- (TMM_S_Flower$mean_CASI_Cg_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO

TMM_S_Flower$TMM_Scrco_Cbp_Co <- (TMM_S_Flower$mean_Cbp_Co_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO # Subgenome Co

TMM_S_Flower$TMM_Scrco_ASI_Co <- (TMM_S_Flower$mean_ASI_Co_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_EUR_Co <- (TMM_S_Flower$mean_EUR_Co_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_ME_Co <- (TMM_S_Flower$mean_ME_Co_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO
TMM_S_Flower$TMM_Scrco_CASI_Co <- (TMM_S_Flower$mean_CASI_Co_F - TMM_S_Flower$Ave_CRCO) / TMM_S_Flower$Ave_CRCO

head(TMM_S_Flower)
# in CG vs CO
TMM_S_Flower$TMM_Scgco_Cbp_Cg <- (TMM_S_Flower$mean_Cbp_Cg_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO  # Subgenome Cg 

TMM_S_Flower$TMM_Scgco_ASI_Cg <- (TMM_S_Flower$mean_ASI_Cg_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_EUR_Cg <- (TMM_S_Flower$mean_EUR_Cg_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_ME_Cg <- (TMM_S_Flower$mean_ME_Cg_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_CASI_Cg <- (TMM_S_Flower$mean_CASI_Cg_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO

TMM_S_Flower$TMM_Scgco_Cbp_Co <- (TMM_S_Flower$mean_Cbp_Co_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO  # Subgenome Co

TMM_S_Flower$TMM_Scgco_ASI_Co <- (TMM_S_Flower$mean_ASI_Co_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_EUR_Co <- (TMM_S_Flower$mean_EUR_Co_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_ME_Co <- (TMM_S_Flower$mean_ME_Co_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO
TMM_S_Flower$TMM_Scgco_CASI_Co <- (TMM_S_Flower$mean_CASI_Co_F - TMM_S_Flower$Ave_CGCO) / TMM_S_Flower$Ave_CGCO

head(TMM_S_Flower)
# in CG vs CR
TMM_S_Flower$TMM_Scgcr_Cbp_Cg <- (TMM_S_Flower$mean_Cbp_Cg_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR  # Subgenome Cg 

TMM_S_Flower$TMM_Scgcr_ASI_Cg <- (TMM_S_Flower$mean_ASI_Cg_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_EUR_Cg <- (TMM_S_Flower$mean_EUR_Cg_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_ME_Cg <- (TMM_S_Flower$mean_ME_Cg_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_CASI_Cg <- (TMM_S_Flower$mean_CASI_Cg_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR

TMM_S_Flower$TMM_Scgcr_Cbp_Co <- (TMM_S_Flower$mean_Cbp_Co_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR  # Subgenome Co

TMM_S_Flower$TMM_Scgcr_ASI_Co <- (TMM_S_Flower$mean_ASI_Co_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_EUR_Co <- (TMM_S_Flower$mean_EUR_Co_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_ME_Co <- (TMM_S_Flower$mean_ME_Co_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR
TMM_S_Flower$TMM_Scgcr_CASI_Co <- (TMM_S_Flower$mean_CASI_Co_F - TMM_S_Flower$Ave_CGCR) / TMM_S_Flower$Ave_CGCR

head(TMM_S_Flower)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_Cbp_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_Cbp_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_EUR_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_EUR_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ME_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ME_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_CASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_CASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_Cbp_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_Cbp_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ASI_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_EUR_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_EUR_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ME_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_ME_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_CASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CR_F > TMM_S_Flower$mean_CO_F),]$TMM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_Cbp_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_Cbp_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_EUR_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_EUR_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ME_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ME_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_CASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_CASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_Cbp_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_Cbp_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ASI_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_EUR_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_EUR_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ME_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_ME_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_CASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CO_F),]$TMM_Scgco_CASI_Co * (-1)



# CG vs CR, CG > CR.

TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_Cbp_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_Cbp_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_EUR_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_EUR_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ME_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ME_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_CASI_Cg <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_CASI_Cg * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_Cbp_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_Cbp_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ASI_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_EUR_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_EUR_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ME_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_ME_Co * (-1)
TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_CASI_Co <- TMM_S_Flower[(TMM_S_Flower$mean_CG_F > TMM_S_Flower$mean_CR_F),]$TMM_Scgcr_CASI_Co * (-1)

# Calculate delta S by comparing the absolute values of TMM_S_Cbp_Cr and TMM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
TMM_S_Flower$delta_TMM_Scrco_Cbp <- abs(TMM_S_Flower$TMM_Scrco_Cbp_Cg) - abs(TMM_S_Flower$TMM_Scrco_Cbp_Co)
TMM_S_Flower$delta_TMM_Scrco_ASI <- abs(TMM_S_Flower$TMM_Scrco_ASI_Cg) - abs(TMM_S_Flower$TMM_Scrco_ASI_Co)
TMM_S_Flower$delta_TMM_Scrco_EUR <- abs(TMM_S_Flower$TMM_Scrco_EUR_Cg) - abs(TMM_S_Flower$TMM_Scrco_EUR_Co)
TMM_S_Flower$delta_TMM_Scrco_ME <- abs(TMM_S_Flower$TMM_Scrco_ME_Cg) - abs(TMM_S_Flower$TMM_Scrco_ME_Co)
TMM_S_Flower$delta_TMM_Scrco_CASI <- abs(TMM_S_Flower$TMM_Scrco_CASI_Cg) - abs(TMM_S_Flower$TMM_Scrco_CASI_Co)
head(TMM_S_Flower)

# CG vs CO
TMM_S_Flower$delta_TMM_Scgco_Cbp <- abs(TMM_S_Flower$TMM_Scgco_Cbp_Cg) - abs(TMM_S_Flower$TMM_Scgco_Cbp_Co)
TMM_S_Flower$delta_TMM_Scgco_ASI <- abs(TMM_S_Flower$TMM_Scgco_ASI_Cg) - abs(TMM_S_Flower$TMM_Scgco_ASI_Co)
TMM_S_Flower$delta_TMM_Scgco_EUR <- abs(TMM_S_Flower$TMM_Scgco_EUR_Cg) - abs(TMM_S_Flower$TMM_Scgco_EUR_Co)
TMM_S_Flower$delta_TMM_Scgco_ME <- abs(TMM_S_Flower$TMM_Scgco_ME_Cg) - abs(TMM_S_Flower$TMM_Scgco_ME_Co)
TMM_S_Flower$delta_TMM_Scgco_CASI <- abs(TMM_S_Flower$TMM_Scgco_CASI_Cg) - abs(TMM_S_Flower$TMM_Scgco_CASI_Co)
head(TMM_S_Flower)

# CG vs CR
TMM_S_Flower$delta_TMM_Scgcr_Cbp <- abs(TMM_S_Flower$TMM_Scgcr_Cbp_Cg) - abs(TMM_S_Flower$TMM_Scgcr_Cbp_Co)
TMM_S_Flower$delta_TMM_Scgcr_ASI <- abs(TMM_S_Flower$TMM_Scgcr_ASI_Cg) - abs(TMM_S_Flower$TMM_Scgcr_ASI_Co)
TMM_S_Flower$delta_TMM_Scgcr_EUR <- abs(TMM_S_Flower$TMM_Scgcr_EUR_Cg) - abs(TMM_S_Flower$TMM_Scgcr_EUR_Co)
TMM_S_Flower$delta_TMM_Scgcr_ME <- abs(TMM_S_Flower$TMM_Scgcr_ME_Cg) - abs(TMM_S_Flower$TMM_Scgcr_ME_Co)
TMM_S_Flower$delta_TMM_Scgcr_CASI <- abs(TMM_S_Flower$TMM_Scgcr_CASI_Cg) - abs(TMM_S_Flower$TMM_Scgcr_CASI_Co)
head(TMM_S_Flower)
# load genes significantly expressed between Cr Co in flower.
DE.F.crco <- read.table("OutputData/DE_CrCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.crco)
dim(DE.F.crco)

DE.F.cgco <- read.table("OutputData/DE_CgCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgco)
dim(DE.F.cgco)

DE.F.cgcr <- read.table("OutputData/DE_CgCr_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgcr)
dim(DE.F.cgcr)

# keep genes with FDR ≤ 0.01
DE.F.crco.keep <- TMM_S_Flower$Gene %in% row.names(DE.F.crco)
TMM_S_Flower$Gene %in% row.names(DE.F.crco)
DE.F.cgco.keep <- TMM_S_Flower$Gene %in% row.names(DE.F.cgco)
TMM_S_Flower$Gene %in% row.names(DE.F.cgco)
# DE.F.cgco.keep <- TMM_S_Flower$Gene %in% row.names(DE.F.crco)
# TMM_S_Flower$Gene %in% row.names(DE.F.crco)
DE.F.cgcr.keep <- TMM_S_Flower$Gene %in% row.names(DE.F.cgcr)
TMM_S_Flower$Gene %in% row.names(DE.F.cgcr)

TMM_Scrco_Flower_FDR <- TMM_S_Flower[DE.F.crco.keep,]
head(TMM_Scrco_Flower_FDR)
dim(TMM_Scrco_Flower_FDR)
TMM_Scgco_Flower_FDR <- TMM_S_Flower[DE.F.cgco.keep,]
head(TMM_Scgco_Flower_FDR)
dim(TMM_Scgco_Flower_FDR)
TMM_Scgcr_Flower_FDR <- TMM_S_Flower[DE.F.cgcr.keep,]
head(TMM_Scgcr_Flower_FDR)
dim(TMM_Scgcr_Flower_FDR)

##################################

# Create a new dataframe for Leaf
TMM_S_Leaf <- data.frame(Gene=rownames(TMM))
TMM_S_Leaf$mean_CG_L <- rowMeans(TMM[,CG_L], na.rm = T) # Average value of CG expression in Leaf
TMM_S_Leaf$mean_CR_L <- rowMeans(TMM[,CR_L], na.rm = T)
TMM_S_Leaf$mean_CO_L <- rowMeans(TMM[,CO_L], na.rm = T)

TMM_S_Leaf$mean_Cbp_Cg_L <- rowMeans(TMM[,Cbp_Cg_L], na.rm = T)
TMM_S_Leaf$mean_Cbp_Co_L <- rowMeans(TMM[,Cbp_Co_L], na.rm = T)

TMM_S_Leaf$mean_ASI_Cg_L <- rowMeans(TMM[,ASI_Cg_L], na.rm = T)
TMM_S_Leaf$mean_EUR_Cg_L <- rowMeans(TMM[,EUR_Cg_L], na.rm = T)
TMM_S_Leaf$mean_ME_Cg_L <- rowMeans(TMM[,ME_Cg_L], na.rm = T)
TMM_S_Leaf$mean_CASI_Cg_L <- rowMeans(TMM[,CASI_Cg_L], na.rm = T)

TMM_S_Leaf$mean_ASI_Co_L <- rowMeans(TMM[,ASI_Co_L], na.rm = T)
TMM_S_Leaf$mean_EUR_Co_L <- rowMeans(TMM[,EUR_Co_L], na.rm = T)
TMM_S_Leaf$mean_ME_Co_L <- rowMeans(TMM[,ME_Co_L], na.rm = T)
TMM_S_Leaf$mean_CASI_Co_L <- rowMeans(TMM[,CASI_Co_L], na.rm = T)

head(TMM[,CG_L])
head(TMM[,EUR_Cg_L])
head(TMM_S_Leaf)

# Average expression level in parents of CR & CO.
TMM_S_Leaf$Ave_CRCO <- rowMeans(TMM_S_Leaf[,c("mean_CR_L","mean_CO_L")])
TMM_S_Leaf$Ave_CGCO <- rowMeans(TMM_S_Leaf[,c("mean_CG_L","mean_CO_L")])
TMM_S_Leaf$Ave_CGCR <- rowMeans(TMM_S_Leaf[,c("mean_CG_L","mean_CR_L")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
TMM_S_Leaf$TMM_Scrco_Cbp_Cg <- (TMM_S_Leaf$mean_Cbp_Cg_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO # Subgenome Cg 

TMM_S_Leaf$TMM_Scrco_ASI_Cg <- (TMM_S_Leaf$mean_ASI_Cg_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_EUR_Cg <- (TMM_S_Leaf$mean_EUR_Cg_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_ME_Cg <- (TMM_S_Leaf$mean_ME_Cg_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_CASI_Cg <- (TMM_S_Leaf$mean_CASI_Cg_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO

TMM_S_Leaf$TMM_Scrco_Cbp_Co <- (TMM_S_Leaf$mean_Cbp_Co_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO # Subgenome Co

TMM_S_Leaf$TMM_Scrco_ASI_Co <- (TMM_S_Leaf$mean_ASI_Co_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_EUR_Co <- (TMM_S_Leaf$mean_EUR_Co_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_ME_Co <- (TMM_S_Leaf$mean_ME_Co_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO
TMM_S_Leaf$TMM_Scrco_CASI_Co <- (TMM_S_Leaf$mean_CASI_Co_L - TMM_S_Leaf$Ave_CRCO) / TMM_S_Leaf$Ave_CRCO

head(TMM_S_Leaf)
# in CG vs CO
TMM_S_Leaf$TMM_Scgco_Cbp_Cg <- (TMM_S_Leaf$mean_Cbp_Cg_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO  # Subgenome Cg 

TMM_S_Leaf$TMM_Scgco_ASI_Cg <- (TMM_S_Leaf$mean_ASI_Cg_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_EUR_Cg <- (TMM_S_Leaf$mean_EUR_Cg_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_ME_Cg <- (TMM_S_Leaf$mean_ME_Cg_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_CASI_Cg <- (TMM_S_Leaf$mean_CASI_Cg_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO

TMM_S_Leaf$TMM_Scgco_Cbp_Co <- (TMM_S_Leaf$mean_Cbp_Co_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO  # Subgenome Co

TMM_S_Leaf$TMM_Scgco_ASI_Co <- (TMM_S_Leaf$mean_ASI_Co_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_EUR_Co <- (TMM_S_Leaf$mean_EUR_Co_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_ME_Co <- (TMM_S_Leaf$mean_ME_Co_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO
TMM_S_Leaf$TMM_Scgco_CASI_Co <- (TMM_S_Leaf$mean_CASI_Co_L - TMM_S_Leaf$Ave_CGCO) / TMM_S_Leaf$Ave_CGCO

head(TMM_S_Leaf)
# in CG vs CR
TMM_S_Leaf$TMM_Scgcr_Cbp_Cg <- (TMM_S_Leaf$mean_Cbp_Cg_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR  # Subgenome Cg 

TMM_S_Leaf$TMM_Scgcr_ASI_Cg <- (TMM_S_Leaf$mean_ASI_Cg_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_EUR_Cg <- (TMM_S_Leaf$mean_EUR_Cg_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_ME_Cg <- (TMM_S_Leaf$mean_ME_Cg_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_CASI_Cg <- (TMM_S_Leaf$mean_CASI_Cg_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR

TMM_S_Leaf$TMM_Scgcr_Cbp_Co <- (TMM_S_Leaf$mean_Cbp_Co_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR  # Subgenome Co

TMM_S_Leaf$TMM_Scgcr_ASI_Co <- (TMM_S_Leaf$mean_ASI_Co_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_EUR_Co <- (TMM_S_Leaf$mean_EUR_Co_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_ME_Co <- (TMM_S_Leaf$mean_ME_Co_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR
TMM_S_Leaf$TMM_Scgcr_CASI_Co <- (TMM_S_Leaf$mean_CASI_Co_L - TMM_S_Leaf$Ave_CGCR) / TMM_S_Leaf$Ave_CGCR

head(TMM_S_Leaf)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_Cbp_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_Cbp_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_EUR_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_EUR_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ME_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ME_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_CASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_CASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_Cbp_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_Cbp_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ASI_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_EUR_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_EUR_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ME_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_ME_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_CASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CR_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scrco_CASI_Co * (-1)

# CG vs CO, CG > CO.

TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_Cbp_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_Cbp_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_EUR_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_EUR_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ME_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ME_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_CASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_CASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_Cbp_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_Cbp_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ASI_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_EUR_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_EUR_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ME_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_ME_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_CASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CO_L),]$TMM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.

TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_Cbp_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_Cbp_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_EUR_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_EUR_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ME_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ME_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_CASI_Cg <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_CASI_Cg * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_Cbp_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_Cbp_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ASI_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_EUR_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_EUR_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ME_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_ME_Co * (-1)
TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_CASI_Co <- TMM_S_Leaf[(TMM_S_Leaf$mean_CG_L > TMM_S_Leaf$mean_CR_L),]$TMM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of TMM_S_Cbp_Cr and TMM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
TMM_S_Leaf$delta_TMM_Scrco_Cbp <- abs(TMM_S_Leaf$TMM_Scrco_Cbp_Co) - abs(TMM_S_Leaf$TMM_Scrco_Cbp_Cg)
TMM_S_Leaf$delta_TMM_Scrco_ASI <- abs(TMM_S_Leaf$TMM_Scrco_ASI_Co) - abs(TMM_S_Leaf$TMM_Scrco_ASI_Cg)
TMM_S_Leaf$delta_TMM_Scrco_EUR <- abs(TMM_S_Leaf$TMM_Scrco_EUR_Co) - abs(TMM_S_Leaf$TMM_Scrco_EUR_Cg)
TMM_S_Leaf$delta_TMM_Scrco_ME <- abs(TMM_S_Leaf$TMM_Scrco_ME_Co) - abs(TMM_S_Leaf$TMM_Scrco_ME_Cg)
TMM_S_Leaf$delta_TMM_Scrco_CASI <- abs(TMM_S_Leaf$TMM_Scrco_CASI_Co) - abs(TMM_S_Leaf$TMM_Scrco_CASI_Cg)
head(TMM_S_Leaf)

# CG vs CO
TMM_S_Leaf$delta_TMM_Scgco_Cbp <- abs(TMM_S_Leaf$TMM_Scgco_Cbp_Co) - abs(TMM_S_Leaf$TMM_Scgco_Cbp_Cg)
TMM_S_Leaf$delta_TMM_Scgco_ASI <- abs(TMM_S_Leaf$TMM_Scgco_ASI_Co) - abs(TMM_S_Leaf$TMM_Scgco_ASI_Cg)
TMM_S_Leaf$delta_TMM_Scgco_EUR <- abs(TMM_S_Leaf$TMM_Scgco_EUR_Co) - abs(TMM_S_Leaf$TMM_Scgco_EUR_Cg)
TMM_S_Leaf$delta_TMM_Scgco_ME <- abs(TMM_S_Leaf$TMM_Scgco_ME_Co) - abs(TMM_S_Leaf$TMM_Scgco_ME_Cg)
TMM_S_Leaf$delta_TMM_Scgco_CASI <- abs(TMM_S_Leaf$TMM_Scgco_CASI_Co) - abs(TMM_S_Leaf$TMM_Scgco_CASI_Cg)
head(TMM_S_Leaf)

# CG vs CR
TMM_S_Leaf$delta_TMM_Scgcr_Cbp <- abs(TMM_S_Leaf$TMM_Scgcr_Cbp_Co) - abs(TMM_S_Leaf$TMM_Scgcr_Cbp_Cg)
TMM_S_Leaf$delta_TMM_Scgcr_ASI <- abs(TMM_S_Leaf$TMM_Scgcr_ASI_Co) - abs(TMM_S_Leaf$TMM_Scgcr_ASI_Cg)
TMM_S_Leaf$delta_TMM_Scgcr_EUR <- abs(TMM_S_Leaf$TMM_Scgcr_EUR_Co) - abs(TMM_S_Leaf$TMM_Scgcr_EUR_Cg)
TMM_S_Leaf$delta_TMM_Scgcr_ME <- abs(TMM_S_Leaf$TMM_Scgcr_ME_Co) - abs(TMM_S_Leaf$TMM_Scgcr_ME_Cg)
TMM_S_Leaf$delta_TMM_Scgcr_CASI <- abs(TMM_S_Leaf$TMM_Scgcr_CASI_Co) - abs(TMM_S_Leaf$TMM_Scgcr_CASI_Cg)
head(TMM_S_Leaf)
# load genes significantly expressed between Cr Co in flower.
DE.L.crco <- read.table("OutputData/DE_CrCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.crco)
dim(DE.L.crco)

DE.L.cgco <- read.table("OutputData/DE_CgCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgco)
dim(DE.L.cgco)

DE.L.cgcr <- read.table("OutputData/DE_CgCr_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgcr)
dim(DE.L.cgcr)

# keep genes with FDR ≤ 0.01
DE.L.crco.keep <- TMM_S_Leaf$Gene %in% row.names(DE.L.crco)
TMM_S_Leaf$Gene %in% row.names(DE.L.crco)
DE.L.cgco.keep <- TMM_S_Leaf$Gene %in% row.names(DE.L.cgco)
# DE.L.cgco.keep <- TMM_S_Leaf$Gene %in% row.names(DE.L.crco)
TMM_S_Leaf$Gene %in% row.names(DE.L.cgco)
DE.L.cgcr.keep <- TMM_S_Leaf$Gene %in% row.names(DE.L.cgcr)
TMM_S_Leaf$Gene %in% row.names(DE.L.cgcr)

TMM_Scrco_Leaf_FDR <- TMM_S_Leaf[DE.L.crco.keep,]
head(TMM_Scrco_Leaf_FDR)
dim(TMM_Scrco_Leaf_FDR)
TMM_Scgco_Leaf_FDR <- TMM_S_Leaf[DE.L.cgco.keep,]
head(TMM_Scgco_Leaf_FDR)
dim(TMM_Scgco_Leaf_FDR)
TMM_Scgcr_Leaf_FDR <- TMM_S_Leaf[DE.L.cgcr.keep,]
head(TMM_Scgcr_Leaf_FDR)
dim(TMM_Scgcr_Leaf_FDR)

#############################

# Create a new dataframe for Root
TMM_S_Root <- data.frame(Gene=rownames(TMM))
TMM_S_Root$mean_CG_R <- rowMeans(TMM[,CG_R], na.rm = T) # Average value of CG expression in Root
TMM_S_Root$mean_CR_R <- rowMeans(TMM[,CR_R], na.rm = T)
TMM_S_Root$mean_CO_R <- rowMeans(TMM[,CO_R], na.rm = T)

TMM_S_Root$mean_Cbp_Cg_R <- rowMeans(TMM[,Cbp_Cg_R], na.rm = T)
TMM_S_Root$mean_Cbp_Co_R <- rowMeans(TMM[,Cbp_Co_R], na.rm = T)

TMM_S_Root$mean_ASI_Cg_R <- rowMeans(TMM[,ASI_Cg_R], na.rm = T)
TMM_S_Root$mean_EUR_Cg_R <- rowMeans(TMM[,EUR_Cg_R], na.rm = T)
TMM_S_Root$mean_ME_Cg_R <- rowMeans(TMM[,ME_Cg_R], na.rm = T)
TMM_S_Root$mean_CASI_Cg_R <- rowMeans(TMM[,CASI_Cg_R], na.rm = T)

TMM_S_Root$mean_ASI_Co_R <- rowMeans(TMM[,ASI_Co_R], na.rm = T)
TMM_S_Root$mean_EUR_Co_R <- rowMeans(TMM[,EUR_Co_R], na.rm = T)
TMM_S_Root$mean_ME_Co_R <- rowMeans(TMM[,ME_Co_R], na.rm = T)
TMM_S_Root$mean_CASI_Co_R <- rowMeans(TMM[,CASI_Co_R], na.rm = T)

head(TMM[,CG_R])
head(TMM[,EUR_Cg_R])
head(TMM_S_Root)

# Average expression level in parents of CR & CO.
TMM_S_Root$Ave_CRCO <- rowMeans(TMM_S_Root[,c("mean_CR_R","mean_CO_R")])
TMM_S_Root$Ave_CGCO <- rowMeans(TMM_S_Root[,c("mean_CG_R","mean_CO_R")])
TMM_S_Root$Ave_CGCR <- rowMeans(TMM_S_Root[,c("mean_CG_R","mean_CR_R")])
# Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
# Where Eij is the average expression of a given transcript i in a given genetic background j 
# CG, CO, CR for parental species, and Cbpcg or Cbpco for subgenomes of C.bursa-pastoris.

# in CR vs CO
TMM_S_Root$TMM_Scrco_Cbp_Cg <- (TMM_S_Root$mean_Cbp_Cg_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO # Subgenome Cg 

TMM_S_Root$TMM_Scrco_ASI_Cg <- (TMM_S_Root$mean_ASI_Cg_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_EUR_Cg <- (TMM_S_Root$mean_EUR_Cg_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_ME_Cg <- (TMM_S_Root$mean_ME_Cg_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_CASI_Cg <- (TMM_S_Root$mean_CASI_Cg_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO

TMM_S_Root$TMM_Scrco_Cbp_Co <- (TMM_S_Root$mean_Cbp_Co_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO # Subgenome Co

TMM_S_Root$TMM_Scrco_ASI_Co <- (TMM_S_Root$mean_ASI_Co_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_EUR_Co <- (TMM_S_Root$mean_EUR_Co_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_ME_Co <- (TMM_S_Root$mean_ME_Co_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO
TMM_S_Root$TMM_Scrco_CASI_Co <- (TMM_S_Root$mean_CASI_Co_R - TMM_S_Root$Ave_CRCO) / TMM_S_Root$Ave_CRCO

head(TMM_S_Root)
# in CG vs CO
TMM_S_Root$TMM_Scgco_Cbp_Cg <- (TMM_S_Root$mean_Cbp_Cg_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO  # Subgenome Cg 

TMM_S_Root$TMM_Scgco_ASI_Cg <- (TMM_S_Root$mean_ASI_Cg_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_EUR_Cg <- (TMM_S_Root$mean_EUR_Cg_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_ME_Cg <- (TMM_S_Root$mean_ME_Cg_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_CASI_Cg <- (TMM_S_Root$mean_CASI_Cg_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO

TMM_S_Root$TMM_Scgco_Cbp_Co <- (TMM_S_Root$mean_Cbp_Co_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO  # Subgenome Co

TMM_S_Root$TMM_Scgco_ASI_Co <- (TMM_S_Root$mean_ASI_Co_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_EUR_Co <- (TMM_S_Root$mean_EUR_Co_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_ME_Co <- (TMM_S_Root$mean_ME_Co_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO
TMM_S_Root$TMM_Scgco_CASI_Co <- (TMM_S_Root$mean_CASI_Co_R - TMM_S_Root$Ave_CGCO) / TMM_S_Root$Ave_CGCO

head(TMM_S_Root)
# in CG vs CR
TMM_S_Root$TMM_Scgcr_Cbp_Cg <- (TMM_S_Root$mean_Cbp_Cg_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR  # Subgenome Cg 

TMM_S_Root$TMM_Scgcr_ASI_Cg <- (TMM_S_Root$mean_ASI_Cg_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_EUR_Cg <- (TMM_S_Root$mean_EUR_Cg_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_ME_Cg <- (TMM_S_Root$mean_ME_Cg_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_CASI_Cg <- (TMM_S_Root$mean_CASI_Cg_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR

TMM_S_Root$TMM_Scgcr_Cbp_Co <- (TMM_S_Root$mean_Cbp_Co_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR  # Subgenome Co

TMM_S_Root$TMM_Scgcr_ASI_Co <- (TMM_S_Root$mean_ASI_Co_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_EUR_Co <- (TMM_S_Root$mean_EUR_Co_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_ME_Co <- (TMM_S_Root$mean_ME_Co_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR
TMM_S_Root$TMM_Scgcr_CASI_Co <- (TMM_S_Root$mean_CASI_Co_R - TMM_S_Root$Ave_CGCR) / TMM_S_Root$Ave_CGCR

head(TMM_S_Root)

# Orientate the S index by expression of parents CR and CO. When CR > CO, S x -1.
# so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
# the expression of that transcript in CR or CO, respectively. 

# CR vs CO, CR > CO.

TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_Cbp_Cg <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_Cbp_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_EUR_Cg <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_EUR_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ME_Cg <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ME_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_CASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_CASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_Cbp_Co <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_Cbp_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ASI_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_EUR_Co <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_EUR_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ME_Co <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_ME_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_CASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CR_R > TMM_S_Root$mean_CO_R),]$TMM_Scrco_CASI_Co * (-1)


# CG vs CO, CG > CO.
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_Cbp_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_Cbp_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_EUR_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_EUR_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ME_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ME_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_CASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_CASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_Cbp_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_Cbp_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ASI_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_EUR_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_EUR_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ME_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_ME_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_CASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CO_R),]$TMM_Scgco_CASI_Co * (-1)


# CG vs CR, CG > CR.
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_Cbp_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_Cbp_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_EUR_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_EUR_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ME_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ME_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_CASI_Cg <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_CASI_Cg * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_Cbp_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_Cbp_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ASI_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_EUR_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_EUR_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ME_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_ME_Co * (-1)
TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_CASI_Co <- TMM_S_Root[(TMM_S_Root$mean_CG_R > TMM_S_Root$mean_CR_R),]$TMM_Scgcr_CASI_Co * (-1)


# Calculate delta S by comparing the absolute values of TMM_S_Cbp_Cr and TMM_S_Cbp_Co.
# (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
# CR vs CO
TMM_S_Root$delta_TMM_Scrco_Cbp <- abs(TMM_S_Root$TMM_Scrco_Cbp_Co) - abs(TMM_S_Root$TMM_Scrco_Cbp_Cg)
TMM_S_Root$delta_TMM_Scrco_ASI <- abs(TMM_S_Root$TMM_Scrco_ASI_Co) - abs(TMM_S_Root$TMM_Scrco_ASI_Cg)
TMM_S_Root$delta_TMM_Scrco_EUR <- abs(TMM_S_Root$TMM_Scrco_EUR_Co) - abs(TMM_S_Root$TMM_Scrco_EUR_Cg)
TMM_S_Root$delta_TMM_Scrco_ME <- abs(TMM_S_Root$TMM_Scrco_ME_Co) - abs(TMM_S_Root$TMM_Scrco_ME_Cg)
TMM_S_Root$delta_TMM_Scrco_CASI <- abs(TMM_S_Root$TMM_Scrco_CASI_Co) - abs(TMM_S_Root$TMM_Scrco_CASI_Cg)
head(TMM_S_Root)

# CG vs CO
TMM_S_Root$delta_TMM_Scgco_Cbp <- abs(TMM_S_Root$TMM_Scgco_Cbp_Co) - abs(TMM_S_Root$TMM_Scgco_Cbp_Cg)
TMM_S_Root$delta_TMM_Scgco_ASI <- abs(TMM_S_Root$TMM_Scgco_ASI_Co) - abs(TMM_S_Root$TMM_Scgco_ASI_Cg)
TMM_S_Root$delta_TMM_Scgco_EUR <- abs(TMM_S_Root$TMM_Scgco_EUR_Co) - abs(TMM_S_Root$TMM_Scgco_EUR_Cg)
TMM_S_Root$delta_TMM_Scgco_ME <- abs(TMM_S_Root$TMM_Scgco_ME_Co) - abs(TMM_S_Root$TMM_Scgco_ME_Cg)
TMM_S_Root$delta_TMM_Scgco_CASI <- abs(TMM_S_Root$TMM_Scgco_CASI_Co) - abs(TMM_S_Root$TMM_Scgco_CASI_Cg)
head(TMM_S_Root)

# CG vs CR
TMM_S_Root$delta_TMM_Scgcr_Cbp <- abs(TMM_S_Root$TMM_Scgcr_Cbp_Co) - abs(TMM_S_Root$TMM_Scgcr_Cbp_Cg)
TMM_S_Root$delta_TMM_Scgcr_ASI <- abs(TMM_S_Root$TMM_Scgcr_ASI_Co) - abs(TMM_S_Root$TMM_Scgcr_ASI_Cg)
TMM_S_Root$delta_TMM_Scgcr_EUR <- abs(TMM_S_Root$TMM_Scgcr_EUR_Co) - abs(TMM_S_Root$TMM_Scgcr_EUR_Cg)
TMM_S_Root$delta_TMM_Scgcr_ME <- abs(TMM_S_Root$TMM_Scgcr_ME_Co) - abs(TMM_S_Root$TMM_Scgcr_ME_Cg)
TMM_S_Root$delta_TMM_Scgcr_CASI <- abs(TMM_S_Root$TMM_Scgcr_CASI_Co) - abs(TMM_S_Root$TMM_Scgcr_CASI_Cg)
head(TMM_S_Root)
# load genes significantly expressed between Cr Co in flower.
DE.R.crco <- read.table("OutputData/DE_CrCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.crco)
dim(DE.R.crco)

DE.R.cgco <- read.table("OutputData/DE_CgCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgco)
dim(DE.R.cgco)

DE.R.cgcr <- read.table("OutputData/DE_CgCr_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgcr)
dim(DE.R.cgcr)

# keep genes with FDR ≤ 0.01
DE.R.crco.keep <- TMM_S_Root$Gene %in% row.names(DE.R.crco)
TMM_S_Root$Gene %in% row.names(DE.R.crco)
DE.R.cgco.keep <- TMM_S_Root$Gene %in% row.names(DE.R.cgco)
# DE.R.cgco.keep <- TMM_S_Root$Gene %in% row.names(DE.R.crco)
TMM_S_Root$Gene %in% row.names(DE.R.cgco)
DE.R.cgcr.keep <- TMM_S_Root$Gene %in% row.names(DE.R.cgcr)
TMM_S_Root$Gene %in% row.names(DE.R.cgcr)

TMM_Scrco_Root_FDR <- TMM_S_Root[DE.R.crco.keep,]
head(TMM_Scrco_Root_FDR)
dim(TMM_Scrco_Root_FDR)
TMM_Scgco_Root_FDR <- TMM_S_Root[DE.R.cgco.keep,]
head(TMM_Scgco_Root_FDR)
dim(TMM_Scgco_Root_FDR)
TMM_Scgcr_Root_FDR <- TMM_S_Root[DE.R.cgcr.keep,]
head(TMM_Scgcr_Root_FDR)
dim(TMM_Scgcr_Root_FDR)
################# Create a new dataframe of S for plot ###################
S_normalisation <- c(rep("CPM", 135), rep("logCPM", 135), rep("TMM", 135))
# S_normalisation <- c(rep("TMM", 45)) # only CGvsCO
S_comparison <- rep(c(rep("CgCo", 45), rep("CrCo", 45), rep("CgCr", 45)), 3) # 3 normalizations
# S_comparison <- c(rep("CgCo", 45))
S_index <- rep(c(rep(c("delta_S","S_Co", "S_Cg"), 45)), 3)
# S_index <- c(rep(c("delta_S","S_Co", "S_Cg"), 15))
# 5 Cbp subgroups (all + asi + eur + me + casi), 3 tissues (F,L,R), 3 comparisons, 3 normalisations
S_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 27)
# S_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 3)
S_tissue <- rep(c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15)), 9)
# S_tissue <- c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15))
# Count median S value for each subpopulation. Median value
###########################################################
CPM_cgco_S_value <- c(abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_Cbp_Co)) - abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_Cbp_Cg)), # Flower all  ##### CPM with cgco 
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_Cbp_Co),
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_Cbp_Cg),
                      abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_ASI_Co)) - abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_ASI_Cg)), # Flower ASI
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_ASI_Co),
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_ASI_Cg),
                      abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_EUR_Co)) - abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_EUR_Cg)), # Flower EUR
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_EUR_Co),
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_EUR_Cg),
                      abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_ME_Co)) - abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_ME_Cg)), # Flower ME
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_ME_Co),
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_ME_Cg),
                      abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_CASI_Co)) - abs(median(CPM_Scgco_Flower_FDR$CPM_Scgco_CASI_Cg)), # Flower CASI
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_CASI_Co),
                      median(CPM_Scgco_Flower_FDR$CPM_Scgco_CASI_Cg),
                      abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_Cbp_Co)) - abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_Cbp_Co),
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_Cbp_Cg),
                      abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ASI_Co)) - abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ASI_Cg)), # Leaf ASI
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ASI_Co),
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ASI_Cg),
                      abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_EUR_Co)) - abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_EUR_Cg)), # Leaf EUR
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_EUR_Co),
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_EUR_Cg),
                      abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ME_Co)) - abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ME_Cg)), # Leaf ME
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ME_Co),
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_ME_Cg),
                      abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_CASI_Co)) - abs(median(CPM_Scgco_Leaf_FDR$CPM_Scgco_CASI_Cg)), # Leaf CASI
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_CASI_Co),
                      median(CPM_Scgco_Leaf_FDR$CPM_Scgco_CASI_Cg),
                      abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_Cbp_Co)) - abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_Cbp_Cg)), # Root all  ##### Root
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_Cbp_Co),
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_Cbp_Cg),
                      abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_ASI_Co)) - abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_ASI_Cg)), # Root ASI
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_ASI_Co),
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_ASI_Cg),
                      abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_EUR_Co)) - abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_EUR_Cg)), # Root EUR
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_EUR_Co),
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_EUR_Cg),
                      abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_ME_Co)) - abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_ME_Cg)), # Root ME
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_ME_Co),
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_ME_Cg),
                      abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_CASI_Co)) - abs(median(CPM_Scgco_Root_FDR$CPM_Scgco_CASI_Cg)), # Root CASI
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_CASI_Co),
                      median(CPM_Scgco_Root_FDR$CPM_Scgco_CASI_Cg)
)

CPM_crco_S_value <- c(abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_Cbp_Co)) - abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_Cbp_Cg)), # Flower all  ##### CPM with crco 
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_Cbp_Co),
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_Cbp_Cg),
                      abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_ASI_Co)) - abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_ASI_Cg)), # Flower ASI
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_ASI_Co),
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_ASI_Cg),
                      abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_EUR_Co)) - abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_EUR_Cg)), # Flower EUR
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_EUR_Co),
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_EUR_Cg),
                      abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_ME_Co)) - abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_ME_Cg)), # Flower ME
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_ME_Co),
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_ME_Cg),
                      abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_CASI_Co)) - abs(median(CPM_Scrco_Flower_FDR$CPM_Scrco_CASI_Cg)), # Flower CASI
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_CASI_Co),
                      median(CPM_Scrco_Flower_FDR$CPM_Scrco_CASI_Cg),
                      abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_Cbp_Co)) - abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_Cbp_Co),
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_Cbp_Cg),
                      abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ASI_Co)) - abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ASI_Cg)), # Leaf ASI
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ASI_Co),
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ASI_Cg),
                      abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_EUR_Co)) - abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_EUR_Cg)), # Leaf EUR
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_EUR_Co),
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_EUR_Cg),
                      abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ME_Co)) - abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ME_Cg)), # Leaf ME
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ME_Co),
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_ME_Cg),
                      abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_CASI_Co)) - abs(median(CPM_Scrco_Leaf_FDR$CPM_Scrco_CASI_Cg)), # Leaf CASI
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_CASI_Co),
                      median(CPM_Scrco_Leaf_FDR$CPM_Scrco_CASI_Cg),
                      abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_Cbp_Co)) - abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_Cbp_Cg)), # Root all  ##### Root
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_Cbp_Co),
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_Cbp_Cg),
                      abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_ASI_Co)) - abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_ASI_Cg)), # Root ASI
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_ASI_Co),
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_ASI_Cg),
                      abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_EUR_Co)) - abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_EUR_Cg)), # Root EUR
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_EUR_Co),
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_EUR_Cg),
                      abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_ME_Co)) - abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_ME_Cg)), # Root ME
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_ME_Co),
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_ME_Cg),
                      abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_CASI_Co)) - abs(median(CPM_Scrco_Root_FDR$CPM_Scrco_CASI_Cg)), # Root CASI
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_CASI_Co),
                      median(CPM_Scrco_Root_FDR$CPM_Scrco_CASI_Cg)
)

CPM_cgcr_S_value <- c(abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_Cbp_Co)) - abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_Cbp_Cg)), # Flower all  ##### CPM with cgcr 
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_Cbp_Co),
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_Cbp_Cg),
                      abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ASI_Co)) - abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ASI_Cg)), # Flower ASI
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ASI_Co),
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ASI_Cg),
                      abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_EUR_Co)) - abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_EUR_Cg)), # Flower EUR
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_EUR_Co),
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_EUR_Cg),
                      abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ME_Co)) - abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ME_Cg)), # Flower ME
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ME_Co),
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_ME_Cg),
                      abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_CASI_Co)) - abs(median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_CASI_Cg)), # Flower CASI
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_CASI_Co),
                      median(CPM_Scgcr_Flower_FDR$CPM_Scgcr_CASI_Cg),
                      abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_Cbp_Co)) - abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_Cbp_Co),
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_Cbp_Cg),
                      abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ASI_Co)) - abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ASI_Cg)), # Leaf ASI
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ASI_Co),
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ASI_Cg),
                      abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_EUR_Co)) - abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_EUR_Cg)), # Leaf EUR
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_EUR_Co),
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_EUR_Cg),
                      abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ME_Co)) - abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ME_Cg)), # Leaf ME
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ME_Co),
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_ME_Cg),
                      abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_CASI_Co)) - abs(median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_CASI_Cg)), # Leaf CASI
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_CASI_Co),
                      median(CPM_Scgcr_Leaf_FDR$CPM_Scgcr_CASI_Cg),
                      abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_Cbp_Co)) - abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_Cbp_Cg)), # Root all  ##### Root
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_Cbp_Co),
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_Cbp_Cg),
                      abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ASI_Co)) - abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ASI_Cg)), # Root ASI
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ASI_Co),
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ASI_Cg),
                      abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_EUR_Co)) - abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_EUR_Cg)), # Root EUR
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_EUR_Co),
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_EUR_Cg),
                      abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ME_Co)) - abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ME_Cg)), # Root ME
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ME_Co),
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_ME_Cg),
                      abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_CASI_Co)) - abs(median(CPM_Scgcr_Root_FDR$CPM_Scgcr_CASI_Cg)), # Root CASI
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_CASI_Co),
                      median(CPM_Scgcr_Root_FDR$CPM_Scgcr_CASI_Cg)
)

logCPM_cgco_S_value <- c(abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_Cbp_Co)) - abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_Cbp_Cg)), # Flower all  ##### logCPM with cgco 
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_Cbp_Co),
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_Cbp_Cg),
                         abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ASI_Co)) - abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ASI_Cg)), # Flower ASI
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ASI_Co),
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ASI_Cg),
                         abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_EUR_Co)) - abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_EUR_Cg)), # Flower EUR
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_EUR_Co),
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_EUR_Cg),
                         abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ME_Co)) - abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ME_Cg)), # Flower ME
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ME_Co),
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_ME_Cg),
                         abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_CASI_Co)) - abs(median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_CASI_Cg)), # Flower CASI
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_CASI_Co),
                         median(logCPM_Scgco_Flower_FDR$logCPM_Scgco_CASI_Cg),
                         abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_Cbp_Co)) - abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_Cbp_Cg)), # Leaf all  ##### Leaf
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_Cbp_Co),
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_Cbp_Cg),
                         abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ASI_Co)) - abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ASI_Cg)), # Leaf ASI
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ASI_Co),
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ASI_Cg),
                         abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_EUR_Co)) - abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_EUR_Cg)), # Leaf EUR
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_EUR_Co),
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_EUR_Cg),
                         abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ME_Co)) - abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ME_Cg)), # Leaf ME
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ME_Co),
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_ME_Cg),
                         abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_CASI_Co)) - abs(median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_CASI_Cg)), # Leaf CASI
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_CASI_Co),
                         median(logCPM_Scgco_Leaf_FDR$logCPM_Scgco_CASI_Cg),
                         abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_Cbp_Co)) - abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_Cbp_Cg)), # Root all  ##### Root
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_Cbp_Co),
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_Cbp_Cg),
                         abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ASI_Co)) - abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ASI_Cg)), # Root ASI
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ASI_Co),
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ASI_Cg),
                         abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_EUR_Co)) - abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_EUR_Cg)), # Root EUR
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_EUR_Co),
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_EUR_Cg),
                         abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ME_Co)) - abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ME_Cg)), # Root ME
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ME_Co),
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_ME_Cg),
                         abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_CASI_Co)) - abs(median(logCPM_Scgco_Root_FDR$logCPM_Scgco_CASI_Cg)), # Root CASI
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_CASI_Co),
                         median(logCPM_Scgco_Root_FDR$logCPM_Scgco_CASI_Cg)
)

logCPM_crco_S_value <- c(abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_Cbp_Co)) - abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_Cbp_Cg)), # Flower all  ##### logCPM with crco 
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_Cbp_Co),
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_Cbp_Cg),
                         abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ASI_Co)) - abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ASI_Cg)), # Flower ASI
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ASI_Co),
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ASI_Cg),
                         abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_EUR_Co)) - abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_EUR_Cg)), # Flower EUR
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_EUR_Co),
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_EUR_Cg),
                         abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ME_Co)) - abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ME_Cg)), # Flower ME
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ME_Co),
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_ME_Cg),
                         abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_CASI_Co)) - abs(median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_CASI_Cg)), # Flower CASI
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_CASI_Co),
                         median(logCPM_Scrco_Flower_FDR$logCPM_Scrco_CASI_Cg),
                         abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_Cbp_Co)) - abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_Cbp_Cg)), # Leaf all  ##### Leaf
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_Cbp_Co),
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_Cbp_Cg),
                         abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ASI_Co)) - abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ASI_Cg)), # Leaf ASI
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ASI_Co),
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ASI_Cg),
                         abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_EUR_Co)) - abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_EUR_Cg)), # Leaf EUR
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_EUR_Co),
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_EUR_Cg),
                         abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ME_Co)) - abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ME_Cg)), # Leaf ME
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ME_Co),
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_ME_Cg),
                         abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_CASI_Co)) - abs(median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_CASI_Cg)), # Leaf CASI
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_CASI_Co),
                         median(logCPM_Scrco_Leaf_FDR$logCPM_Scrco_CASI_Cg),
                         abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_Cbp_Co)) - abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_Cbp_Cg)), # Root all  ##### Root
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_Cbp_Co),
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_Cbp_Cg),
                         abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ASI_Co)) - abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ASI_Cg)), # Root ASI
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ASI_Co),
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ASI_Cg),
                         abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_EUR_Co)) - abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_EUR_Cg)), # Root EUR
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_EUR_Co),
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_EUR_Cg),
                         abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ME_Co)) - abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ME_Cg)), # Root ME
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ME_Co),
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_ME_Cg),
                         abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_CASI_Co)) - abs(median(logCPM_Scrco_Root_FDR$logCPM_Scrco_CASI_Cg)), # Root CASI
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_CASI_Co),
                         median(logCPM_Scrco_Root_FDR$logCPM_Scrco_CASI_Cg)
)

logCPM_cgcr_S_value <- c(abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_Cbp_Co)) - abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_Cbp_Cg)), # Flower all  ##### logCPM with cgcr 
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_Cbp_Co),
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_Cbp_Cg),
                         abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ASI_Co)) - abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ASI_Cg)), # Flower ASI
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ASI_Co),
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ASI_Cg),
                         abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_EUR_Co)) - abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_EUR_Cg)), # Flower EUR
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_EUR_Co),
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_EUR_Cg),
                         abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ME_Co)) - abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ME_Cg)), # Flower ME
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ME_Co),
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_ME_Cg),
                         abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_CASI_Co)) - abs(median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_CASI_Cg)), # Flower CASI
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_CASI_Co),
                         median(logCPM_Scgcr_Flower_FDR$logCPM_Scgcr_CASI_Cg),
                         abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_Cbp_Co)) - abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_Cbp_Cg)), # Leaf all  ##### Leaf
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_Cbp_Co),
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_Cbp_Cg),
                         abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ASI_Co)) - abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ASI_Cg)), # Leaf ASI
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ASI_Co),
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ASI_Cg),
                         abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_EUR_Co)) - abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_EUR_Cg)), # Leaf EUR
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_EUR_Co),
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_EUR_Cg),
                         abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ME_Co)) - abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ME_Cg)), # Leaf ME
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ME_Co),
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_ME_Cg),
                         abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_CASI_Co)) - abs(median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_CASI_Cg)), # Leaf CASI
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_CASI_Co),
                         median(logCPM_Scgcr_Leaf_FDR$logCPM_Scgcr_CASI_Cg),
                         abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_Cbp_Co)) - abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_Cbp_Cg)), # Root all  ##### Root
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_Cbp_Co),
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_Cbp_Cg),
                         abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ASI_Co)) - abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ASI_Cg)), # Root ASI
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ASI_Co),
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ASI_Cg),
                         abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_EUR_Co)) - abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_EUR_Cg)), # Root EUR
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_EUR_Co),
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_EUR_Cg),
                         abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ME_Co)) - abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ME_Cg)), # Root ME
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ME_Co),
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_ME_Cg),
                         abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_CASI_Co)) - abs(median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_CASI_Cg)), # Root CASI
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_CASI_Co),
                         median(logCPM_Scgcr_Root_FDR$logCPM_Scgcr_CASI_Cg)
)

TMM_cgco_S_value <- c(abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_Cbp_Co)) - abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_Cbp_Cg)), # Flower all  ##### TMM with cgco 
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_Cbp_Co),
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_Cbp_Cg),
                      abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_ASI_Co)) - abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_ASI_Cg)), # Flower ASI
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_ASI_Co),
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_ASI_Cg),
                      abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_EUR_Co)) - abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_EUR_Cg)), # Flower EUR
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_EUR_Co),
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_EUR_Cg),
                      abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_ME_Co)) - abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_ME_Cg)), # Flower ME
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_ME_Co),
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_ME_Cg),
                      abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_CASI_Co)) - abs(median(TMM_Scgco_Flower_FDR$TMM_Scgco_CASI_Cg)), # Flower CASI
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_CASI_Co),
                      median(TMM_Scgco_Flower_FDR$TMM_Scgco_CASI_Cg),
                      abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_Cbp_Co)) - abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_Cbp_Co),
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_Cbp_Cg),
                      abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ASI_Co)) - abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ASI_Cg)), # Leaf ASI
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ASI_Co),
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ASI_Cg),
                      abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_EUR_Co)) - abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_EUR_Cg)), # Leaf EUR
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_EUR_Co),
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_EUR_Cg),
                      abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ME_Co)) - abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ME_Cg)), # Leaf ME
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ME_Co),
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_ME_Cg),
                      abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_CASI_Co)) - abs(median(TMM_Scgco_Leaf_FDR$TMM_Scgco_CASI_Cg)), # Leaf CASI
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_CASI_Co),
                      median(TMM_Scgco_Leaf_FDR$TMM_Scgco_CASI_Cg),
                      abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_Cbp_Co)) - abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_Cbp_Cg)), # Root all  ##### Root
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_Cbp_Co),
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_Cbp_Cg),
                      abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_ASI_Co)) - abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_ASI_Cg)), # Root ASI
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_ASI_Co),
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_ASI_Cg),
                      abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_EUR_Co)) - abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_EUR_Cg)), # Root EUR
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_EUR_Co),
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_EUR_Cg),
                      abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_ME_Co)) - abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_ME_Cg)), # Root ME
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_ME_Co),
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_ME_Cg),
                      abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_CASI_Co)) - abs(median(TMM_Scgco_Root_FDR$TMM_Scgco_CASI_Cg)), # Root CASI
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_CASI_Co),
                      median(TMM_Scgco_Root_FDR$TMM_Scgco_CASI_Cg)
)

TMM_crco_S_value <- c(abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_Cbp_Co)) - abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_Cbp_Cg)), # Flower all  ##### TMM with crco 
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_Cbp_Co),
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_Cbp_Cg),
                      abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_ASI_Co)) - abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_ASI_Cg)), # Flower ASI
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_ASI_Co),
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_ASI_Cg),
                      abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_EUR_Co)) - abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_EUR_Cg)), # Flower EUR
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_EUR_Co),
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_EUR_Cg),
                      abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_ME_Co)) - abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_ME_Cg)), # Flower ME
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_ME_Co),
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_ME_Cg),
                      abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_CASI_Co)) - abs(median(TMM_Scrco_Flower_FDR$TMM_Scrco_CASI_Cg)), # Flower CASI
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_CASI_Co),
                      median(TMM_Scrco_Flower_FDR$TMM_Scrco_CASI_Cg),
                      abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_Cbp_Co)) - abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_Cbp_Co),
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_Cbp_Cg),
                      abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ASI_Co)) - abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ASI_Cg)), # Leaf ASI
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ASI_Co),
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ASI_Cg),
                      abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_EUR_Co)) - abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_EUR_Cg)), # Leaf EUR
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_EUR_Co),
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_EUR_Cg),
                      abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ME_Co)) - abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ME_Cg)), # Leaf ME
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ME_Co),
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_ME_Cg),
                      abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_CASI_Co)) - abs(median(TMM_Scrco_Leaf_FDR$TMM_Scrco_CASI_Cg)), # Leaf CASI
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_CASI_Co),
                      median(TMM_Scrco_Leaf_FDR$TMM_Scrco_CASI_Cg),
                      abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_Cbp_Co)) - abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_Cbp_Cg)), # Root all  ##### Root
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_Cbp_Co),
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_Cbp_Cg),
                      abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_ASI_Co)) - abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_ASI_Cg)), # Root ASI
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_ASI_Co),
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_ASI_Cg),
                      abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_EUR_Co)) - abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_EUR_Cg)), # Root EUR
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_EUR_Co),
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_EUR_Cg),
                      abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_ME_Co)) - abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_ME_Cg)), # Root ME
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_ME_Co),
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_ME_Cg),
                      abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_CASI_Co)) - abs(median(TMM_Scrco_Root_FDR$TMM_Scrco_CASI_Cg)), # Root CASI
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_CASI_Co),
                      median(TMM_Scrco_Root_FDR$TMM_Scrco_CASI_Cg)
)

TMM_cgcr_S_value <- c(abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_Cbp_Co)) - abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_Cbp_Cg)), # Flower all  ##### TMM with cgcr 
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_Cbp_Co),
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_Cbp_Cg),
                      abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ASI_Co)) - abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ASI_Cg)), # Flower ASI
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ASI_Co),
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ASI_Cg),
                      abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_EUR_Co)) - abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_EUR_Cg)), # Flower EUR
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_EUR_Co),
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_EUR_Cg),
                      abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ME_Co)) - abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ME_Cg)), # Flower ME
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ME_Co),
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_ME_Cg),
                      abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_CASI_Co)) - abs(median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_CASI_Cg)), # Flower CASI
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_CASI_Co),
                      median(TMM_Scgcr_Flower_FDR$TMM_Scgcr_CASI_Cg),
                      abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_Cbp_Co)) - abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_Cbp_Cg)), # Leaf all  ##### Leaf
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_Cbp_Co),
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_Cbp_Cg),
                      abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ASI_Co)) - abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ASI_Cg)), # Leaf ASI
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ASI_Co),
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ASI_Cg),
                      abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_EUR_Co)) - abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_EUR_Cg)), # Leaf EUR
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_EUR_Co),
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_EUR_Cg),
                      abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ME_Co)) - abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ME_Cg)), # Leaf ME
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ME_Co),
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_ME_Cg),
                      abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_CASI_Co)) - abs(median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_CASI_Cg)), # Leaf CASI
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_CASI_Co),
                      median(TMM_Scgcr_Leaf_FDR$TMM_Scgcr_CASI_Cg),
                      abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_Cbp_Co)) - abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_Cbp_Cg)), # Root all  ##### Root
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_Cbp_Co),
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_Cbp_Cg),
                      abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ASI_Co)) - abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ASI_Cg)), # Root ASI
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ASI_Co),
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ASI_Cg),
                      abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_EUR_Co)) - abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_EUR_Cg)), # Root EUR
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_EUR_Co),
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_EUR_Cg),
                      abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ME_Co)) - abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ME_Cg)), # Root ME
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ME_Co),
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_ME_Cg),
                      abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_CASI_Co)) - abs(median(TMM_Scgcr_Root_FDR$TMM_Scgcr_CASI_Cg)), # Root CASI
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_CASI_Co),
                      median(TMM_Scgcr_Root_FDR$TMM_Scgcr_CASI_Cg)
)
###########################################################

S_value <- c(CPM_cgco_S_value, CPM_crco_S_value, CPM_cgcr_S_value,
             logCPM_cgco_S_value, logCPM_crco_S_value, logCPM_cgcr_S_value,
             TMM_cgco_S_value, TMM_crco_S_value, TMM_cgcr_S_value)
# S_value <- TMM_cgco_S_value
# Combine all vectors to a dataframe
S_plot <- cbind(S_normalisation, S_comparison, S_index, S_population, S_tissue, S_value)
# change the format of the dataframe
S_plot <- print.data.frame(data.frame(S_plot), quote=FALSE)
str(S_plot)
S_plot$S_normalisation <- factor(S_plot$S_normalisation, levels = c("CPM", "logCPM", "TMM"))
S_plot$S_comparison <- factor(S_plot$S_comparison, levels = c("CgCo", "CrCo", "CgCr"))
S_plot$S_index <- factor(S_plot$S_index, levels = c("S_Co","S_Cg","delta_S"))
S_plot$S_population <- factor(S_plot$S_population, levels = c("Cbp","Cbp_ASI","Cbp_EUR","Cbp_ME","Cbp_CASI"))
S_plot$S_tissue <- factor(S_plot$S_tissue, levels = c("Flower","Leaf","Root"))
S_plot$S_value <- as.numeric(S_plot$S_value)
str(S_plot)

# Plot
library(ggplot2)
ggplot(S_plot, aes(x=S_tissue, y=S_value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=S_index, shape=S_index),size=2) +
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  #scale_y_continuous(limits=c(-0.15, 0.15)) +
  facet_grid(S_comparison + S_normalisation ~ S_population) + 
  labs(x ="Tissues", y = "Median Si") +
  theme_bw()
#theme(panel.background = element_rect(fill = "white", colour = "grey50"))

library(ggpubr)
ggscatter(S_plot, x = "tissue", y = "value", color = "index", shape = "index",
          ylim = c(-0.3,0.3))


