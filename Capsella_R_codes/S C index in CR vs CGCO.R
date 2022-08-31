###############################################################
#      S and C index of CR in comparison with CG CO           #
###############################################################

# To quantify the similarity between each subgenome expression level and the expression level in the parental 
# species, we developed a similarity index (S). For each transcript i and each subgenome j in {CbpCg, CbpCo}, 
# S was computed as the subgenome relative expression deviation from the mean expression level in the parental
# species, u(i) = (E(ico) + E(icg)) / 2:
# S(ij) = (E(ij) - u(i)) / u(i)

## Set work directory and load files
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

library(ggplot2)

tissue <- c("F", "L", "R")

for(i in tissue){
  print(i)
  
  # set a sub loop for viable name
  if (i == "F"){
    Tissue <- "Flower"
  }else{
    if (i == "L"){
      Tissue <- "Leaf"
    }else{
      Tissue <- "Root"
    }
  }
  
  TMM <- read.table(paste0("OutputData/TMM_CGCRCO_", i, ".txt"), header = T, sep = "\t")
  TMM[TMM == 0] <- NA # substitute na to 0 for the next mean calculation
  
  # Similarity index
  # Create a new dataframe for Flower
  S_CR <- data.frame(Gene=rownames(TMM))
  S_CR$mean_CG <- rowMeans(TMM[,which(sub("[0-9]_.*", "", colnames(TMM)) == "CG")], na.rm = T) # Average value of CG expression in Flower
  S_CR$mean_CR <- rowMeans(TMM[,which(sub("[0-9]_.*", "", colnames(TMM)) == "CR")], na.rm = T)
  S_CR$mean_CO <- rowMeans(TMM[,which(sub("[0-9]_.*", "", colnames(TMM)) == "CO")], na.rm = T)
  
  # Average expression level in parents of CG & CO.
  S_CR$Ave_CGCO <- rowMeans(S_CR[,c("mean_CG","mean_CO")])
  
  # Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
  # Where Eij is the average expression of a given transcript i in a given genetic background j 
  # CG, CO for parental species, and CR as offsprings.
  S_CR$S_index <- (S_CR$mean_CR - S_CR$Ave_CGCO) / S_CR$Ave_CGCO
  
  # Orientate the S index by expression of parents CG and CO. When CG > CO, S x -1.
  # so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
  # the expression of that transcript in CG or CO, respectively. 
  orientate <- S_CR$mean_CG > S_CR$mean_CO
  S_CR[orientate,]$S_index <- S_CR[orientate,]$S_index * (-1)
  
  # count proportion of S < 0
  prop <- length(which(S_CR$S_index < 0)) / length(S_CR$S_index)
  # binormal test for S<0, S>0
  b_test<-binom.test(c(length(which(S_CR$S_index < 0)),length(which(S_CR$S_index > 0))),p=0.5)
  if(b_test$p.value>=0.05){ # set the left label
    leftOut<-paste(round(prop,2)*100,"% ns",sep=" ")}else{
    if(b_test$p.value<0.001){leftOut<-paste(round(prop,2)*100,"% ***",sep="")}else{
      if(b_test$p.value<0.01){leftOut<-paste(round(prop,2)*100,"% **",sep="")}else{
        leftOut<-paste(round(prop,2),"*",sep=" ")}}
    }
  rightOut <- paste0(round(1-prop,2)*100, "%") # set the right label
  myMedian <- paste0("median=",round(median(S_CR$S_index),2))
  
  # plot distribution of S
  p1 <- ggplot(S_CR, aes(x=S_index)) + 
    geom_histogram(aes(y=..density..),
                   color="black", fill="white", bins = 60) +
    geom_vline(aes(xintercept=0),
               color="gray", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=median(S_index)),
               color="red", linetype="dashed", size=1) +
    xlim(c(-1.5, 1.5)) + ylim(c(0,2)) +
    labs(x =expression(S[Co]), y = "Density", title = paste0(Tissue, " All (n=", dim(TMM)[1], ")")) +
    annotate(geom="text", x=-1, y=2, label="Towards CG") +
    annotate(geom="text", x=1, y=2, label="Towards CO") +
    annotate(geom="text", x=-1, y=2*0.9, label=leftOut) +
    annotate(geom="text", x=1, y=2*0.9, label=rightOut) +
    annotate(geom="text", x=-1, y=2*0.8, label=myMedian, color="red") +
    theme_classic()
  assign(paste0("p_", i, "_all"), p1)
  
  # DE genes in different comparisons
  for (com in c("CGCO", "CGCR", "CRCO")) {
    # load DE genes
    DE <- read.table(paste0("OutputData/DE_", com, "_", i, "_FC2FDR.05.txt"), header = T, sep = "\t")
    # keep DE genes
    DE.keep <- S_CR$Gene %in% row.names(DE)
    S_CR_DE <- S_CR[DE.keep,]
    
    # count proportion of S < 0
    prop <- length(which(S_CR_DE$S_index < 0)) / length(S_CR_DE$S_index)
    # binormal test for S<0, S>0
    b_test<-binom.test(c(length(which(S_CR_DE$S_index < 0)),length(which(S_CR_DE$S_index > 0))),p=0.5) #p: hypothesized probability of success
    if(b_test$p.value>=0.05){ # set the left label
      leftOut<-paste(round(prop,2)*100,"% ns",sep=" ")}else{
        if(b_test$p.value<0.001){leftOut<-paste(round(prop,2)*100,"% ***",sep="")}else{
          if(b_test$p.value<0.01){leftOut<-paste(round(prop,2)*100,"% **",sep="")}else{
            leftOut<-paste(round(prop,2),"*",sep=" ")}}
      }
    rightOut <- paste0(round(1-prop,2)*100, "%") # set the right label
    myMedian <- paste0("median=",round(median(S_CR_DE$S_index),2))
    
    if((i == "L") & (com == "CGCR")){
      message(i, "_", com)
      myylim <- c(0,3)
      myPosition <- 3
    }else{
      myylim <- c(0,2)
      myPosition <- 2
    }
    
    if(median(S_CR_DE$S_index) < 0){
      p2 <- ggplot(S_CR_DE, aes(x=S_index)) + 
        geom_histogram(aes(y=..density..),
                       color="black", fill="white", bins = 60) +
        geom_vline(aes(xintercept=0),
                   color="gray", linetype="dashed", size=1) +
        geom_vline(aes(xintercept=median(S_index)),
                   color="red", linetype="dashed", size=1) +
        xlim(c(-1.5, 1.5)) + ylim(myylim) +
        labs(x =expression(S[CR]), y = "Density", title = paste0(Tissue, " DEs of ", com, " (n=", dim(S_CR_DE)[1], ")")) +
        annotate(geom="text", x=-1, y=myPosition, label="Towards CG") +
        annotate(geom="text", x=1, y=myPosition, label="Towards CO") +
        annotate(geom="text", x=-1, y=myPosition*0.9, label=leftOut) +
        annotate(geom="text", x=1, y=myPosition*0.9, label=rightOut) +
        annotate(geom="text", x=-1, y=myPosition*0.8, label=myMedian, color="red") +
        theme_classic()
      assign(paste0("p_", i, "_de_", com), p2)
    }else{
      p2 <- ggplot(S_CR_DE, aes(x=S_index)) + 
        geom_histogram(aes(y=..density..),
                       color="black", fill="white", bins = 60) +
        geom_vline(aes(xintercept=0),
                   color="gray", linetype="dashed", size=1) +
        geom_vline(aes(xintercept=median(S_index)),
                   color="blue", linetype="dashed", size=1) +
        xlim(c(-1.5, 1.5)) + ylim(myylim) +
        labs(x =expression(S[CR]), y = "Density", title = paste0(Tissue, " DEs of ", com, " (n=", dim(S_CR_DE)[1], ")")) +
        annotate(geom="text", x=-1, y=myPosition, label="Towards CG") +
        annotate(geom="text", x=1, y=myPosition, label="Towards CO") +
        annotate(geom="text", x=-1, y=myPosition*0.9, label=leftOut) +
        annotate(geom="text", x=1, y=myPosition*0.9, label=rightOut) +
        annotate(geom="text", x=1, y=myPosition*0.8, label=myMedian, color="blue") +
        theme_classic()
      assign(paste0("p_", i, "_de_", com), p2)
    }
    
  }
}

library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_all, vp = vplayout(1,1))
print(p_L_all, vp = vplayout(1,2))
print(p_R_all, vp = vplayout(1,3))
print(p_F_de_CRCO, vp = vplayout(2,1))
print(p_L_de_CRCO, vp = vplayout(2,2))
print(p_R_de_CRCO, vp = vplayout(2,3))
print(p_F_de_CGCR, vp = vplayout(3,1))
print(p_L_de_CGCR, vp = vplayout(3,2))
print(p_R_de_CGCR, vp = vplayout(3,3))
print(p_F_de_CGCO, vp = vplayout(4,1))
print(p_L_de_CGCO, vp = vplayout(4,2))
print(p_R_de_CGCO, vp = vplayout(4,3))



# keep genes with FDR ≤ 0.01
DE.F.cgco.keep <- S_CR$Gene %in% row.names(DE.F.cgco)
S_CR$Gene %in% row.names(DE.F.cgco)

S_CR_F.cgco <- S_CR[DE.F.cgco.keep,]
head(S_CR_F.cgco)
dim(S_CR_F.cgco)

p4 <- ggplot(S_CR_F.cgco, aes(x=S_CR_F)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_F)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Flower DE_CGCO") +
  theme_classic()

## DE.L.cgco
DE.L.cgco <- read.table("OutputData/DE_CgCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgco)
dim(DE.L.cgco)

# keep genes with FDR ≤ 0.01
DE.L.cgco.keep <- S_CR$Gene %in% row.names(DE.L.cgco)
S_CR$Gene %in% row.names(DE.L.cgco)

S_CR_L.cgco <- S_CR[DE.L.cgco.keep,]
head(S_CR_L.cgco)
dim(S_CR_L.cgco)

p5 <- ggplot(S_CR_L.cgco, aes(x=S_CR_L)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_L)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "leaf DE_CGCO") +
  theme_classic()

## DE.R.cgco
DE.R.cgco <- read.table("OutputData/DE_CgCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgco)
dim(DE.R.cgco)

# keep genes with FDR ≤ 0.01
DE.R.cgco.keep <- S_CR$Gene %in% row.names(DE.R.cgco)
S_CR$Gene %in% row.names(DE.R.cgco)

S_CR_R.cgco <- S_CR[DE.R.cgco.keep,]
head(S_CR_R.cgco)
dim(S_CR_R.cgco)

p6 <- ggplot(S_CR_R.cgco, aes(x=S_CR_R)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_R)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Root DE_CGCO") +
  theme_classic()

# load genes significantly expressed between CrCo in flower.
DE.F.crco <- read.table("OutputData/DE_CrCo_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.crco)
dim(DE.F.crco)

# keep genes with FDR ≤ 0.01
DE.F.crco.keep <- S_CR$Gene %in% row.names(DE.F.crco)
S_CR$Gene %in% row.names(DE.F.crco)

S_CR_F.crco <- S_CR[DE.F.crco.keep,]
head(S_CR_F.crco)
dim(S_CR_F.crco)

p7 <- ggplot(S_CR_F.crco, aes(x=S_CR_F)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_F)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Flower DE_CRCO") +
  theme_classic()

## DE.L.crco
DE.L.crco <- read.table("OutputData/DE_CrCo_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.crco)
dim(DE.L.crco)

# keep genes with FDR ≤ 0.01
DE.L.crco.keep <- S_CR$Gene %in% row.names(DE.L.crco)
S_CR$Gene %in% row.names(DE.L.crco)

S_CR_L.crco <- S_CR[DE.L.crco.keep,]
head(S_CR_L.crco)
dim(S_CR_L.crco)

p8 <- ggplot(S_CR_L.crco, aes(x=S_CR_L)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_L)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "leaf DE_CRCO") +
  theme_classic()

## DE.R.crco
DE.R.crco <- read.table("OutputData/DE_CrCo_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.crco)
dim(DE.R.crco)

# keep genes with FDR ≤ 0.01
DE.R.crco.keep <- S_CR$Gene %in% row.names(DE.R.crco)
S_CR$Gene %in% row.names(DE.R.crco)

S_CR_R.crco <- S_CR[DE.R.crco.keep,]
head(S_CR_R.crco)
dim(S_CR_R.crco)

p9 <- ggplot(S_CR_R.crco, aes(x=S_CR_R)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_R)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Root DE_CRCO") +
  theme_classic()

# load genes significantly expressed between CgCr in flower.
DE.F.cgcr <- read.table("OutputData/DE_CgCr_F genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.F.cgcr)
dim(DE.F.cgcr)

# keep genes with FDR ≤ 0.01
DE.F.cgcr.keep <- S_CR$Gene %in% row.names(DE.F.cgcr)
S_CR$Gene %in% row.names(DE.F.cgcr)

S_CR_F.cgcr <- S_CR[DE.F.cgcr.keep,]
head(S_CR_F.cgcr)
dim(S_CR_F.cgcr)

p10 <- ggplot(S_CR_F.cgcr, aes(x=S_CR_F)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_F)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Flower DE_CGCR") +
  theme_classic()

## DE.L.cgcr
DE.L.cgcr <- read.table("OutputData/DE_CgCr_L genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.L.cgcr)
dim(DE.L.cgcr)

# keep genes with FDR ≤ 0.01
DE.L.cgcr.keep <- S_CR$Gene %in% row.names(DE.L.cgcr)
S_CR$Gene %in% row.names(DE.L.cgcr)

S_CR_L.cgcr <- S_CR[DE.L.cgcr.keep,]
head(S_CR_L.cgcr)
dim(S_CR_L.cgcr)

p11<- ggplot(S_CR_L.cgcr, aes(x=S_CR_L)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_L)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "leaf DE_CGCR") +
  theme_classic()

## DE.R.cgcr
DE.R.cgcr <- read.table("OutputData/DE_CgCr_R genes FDR 0.01.txt", header = T, sep = "\t")
head(DE.R.cgcr)
dim(DE.R.cgcr)

# keep genes with FDR ≤ 0.01
DE.R.cgcr.keep <- S_CR$Gene %in% row.names(DE.R.cgcr)
S_CR$Gene %in% row.names(DE.R.cgcr)

S_CR_R.cgcr <- S_CR[DE.R.cgcr.keep,]
head(S_CR_R.cgcr)
dim(S_CR_R.cgcr)

p12 <- ggplot(S_CR_R.cgcr, aes(x=S_CR_R)) + 
  geom_histogram(aes(y=..density..),
                 color="black", fill="white", bins = 60) +
  geom_vline(aes(xintercept=0),
             color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(S_CR_R)),
             color="blue", linetype="dashed", size=1) +
  xlim(c(-1, 1)) +
  labs(x ="S CR", y = "Density", title = "Root DE_CGCR") +
  theme_classic()
library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p3, vp = vplayout(1,3))
print(p4, vp = vplayout(2,1))
print(p5, vp = vplayout(2,2))
print(p6, vp = vplayout(2,3))
print(p11, vp = vplayout(2,1))
print(p12, vp = vplayout(2,2))
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
write.table(Scrco_Flower_FDR, file = "OutputData/S and C indices in Flower FDR genes.txt", quote = F, sep = "\t", row.names = F)
write.table(Scrco_Leaf_FDR, file = "OutputData/S and C indices in Leaf FDR genes.txt", quote = F, sep = "\t", row.names = F)
write.table(Scrco_Root_FDR, file = "OutputData/S and C indices in Root FDR genes.txt", quote = F, sep = "\t", row.names = F)

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




