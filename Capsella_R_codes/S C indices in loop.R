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

TMM <- read.table("InputData/Capsella_averageTMM_FLR.txt", header = T, sep = "\t")
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
########################################################################################################
# Similarity index -- CG vs CO
tissue<-c("F","L","R") # set tissue

# Loop start
for (i in tissue) {
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
  
  # subData by tissue
  reads <- TMM[, which(sub(".*_","",colnames(TMM)) == i |
                         sub(".*[0-9]_","",colnames(TMM)) == paste(i,"_Cg",sep="") |
                         sub(".*[0-9]_","",colnames(TMM)) == paste(i,"_Co",sep="") )]
  
  # Create a new dataframe for sub-data
  S <- data.frame(Gene=rownames(reads))
  
  S$mean_CG <- rowMeans(reads[,which(sub("[0-9]_.*", "", colnames(reads)) == "CG")], na.rm = T) # Average value of CG expression in sub-data
  S$mean_CR <- rowMeans(reads[,which(sub("[0-9]_.*", "", colnames(reads)) == "CR")], na.rm = T)
  S$mean_CO <- rowMeans(reads[,which(sub("[0-9]_.*", "", colnames(reads)) == "CO")], na.rm = T)
  
  S$mean_Cbp_Cg <- rowMeans(reads[,which(sub(".*_", "", colnames(reads)) == "Cg")], na.rm = T)
  S$mean_Cbp_Co <- rowMeans(reads[,which(sub(".*_", "", colnames(reads)) == "Co")], na.rm = T)
  
  S$mean_ASI_Cg <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "ASI_Cg")], na.rm = T)
  S$mean_EUR_Cg <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "EUR_Cg")], na.rm = T)
  S$mean_ME_Cg <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "ME_Cg")], na.rm = T)
  S$mean_CASI_Cg <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "CASI_Cg")], na.rm = T)
  
  S$mean_ASI_Co <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "ASI_Co")], na.rm = T)
  S$mean_EUR_Co <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "EUR_Co")], na.rm = T)
  S$mean_ME_Co <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "ME_Co")], na.rm = T)
  S$mean_CASI_Co <- rowMeans(reads[,which(sub("[0-9]_.", "", colnames(reads)) == "CASI_Co")], na.rm = T)
  
  # Average expression level in parents of CG & CO.
  S$Ave_CGCO <- rowMeans(S[,c("mean_CG","mean_CO")])
  # Calculate Similarity index: S = (Eij - Ave_CRCO) / Ave_CRCO
  # Where Eij is the average expression of a given transcript i in a given genetic background j 
  # CG, CO for parental species, and CR as offsprings.
  S$Cbp_Cg <- (S$mean_Cbp_Cg - S$Ave_CGCO) / S$Ave_CGCO  # Subgenome Cg 
  S$ASI_Cg <- (S$mean_ASI_Cg - S$Ave_CGCO) / S$Ave_CGCO
  S$EUR_Cg <- (S$mean_EUR_Cg - S$Ave_CGCO) / S$Ave_CGCO
  S$ME_Cg <- (S$mean_ME_Cg - S$Ave_CGCO) / S$Ave_CGCO
  S$CASI_Cg <- (S$mean_CASI_Cg - S$Ave_CGCO) / S$Ave_CGCO
  
  S$Cbp_Co <- (S$mean_Cbp_Co - S$Ave_CGCO) / S$Ave_CGCO  # Subgenome Co
  S$ASI_Co <- (S$mean_ASI_Co - S$Ave_CGCO) / S$Ave_CGCO
  S$EUR_Co <- (S$mean_EUR_Co - S$Ave_CGCO) / S$Ave_CGCO
  S$ME_Co <- (S$mean_ME_Co - S$Ave_CGCO) / S$Ave_CGCO
  S$CASI_Co <- (S$mean_CASI_Co - S$Ave_CGCO) / S$Ave_CGCO
  
  # Orientate the S index by expression of parents CG and CO. When CG > CO, S x -1.
  # so that if Sij < 0 or Sij > 0, the expression of a given transcript in a given subgenome is more similar to 
  # the expression of that transcript in CG or CO, respectively. 
  orientate <- S$mean_CG > S$mean_CO
  
  S[orientate,]$Cbp_Cg <- S[orientate,]$Cbp_Cg * (-1)
  S[orientate,]$ASI_Cg <- S[orientate,]$ASI_Cg * (-1)
  S[orientate,]$EUR_Cg <- S[orientate,]$EUR_Cg * (-1)
  S[orientate,]$ME_Cg <- S[orientate,]$ME_Cg * (-1)
  S[orientate,]$CASI_Cg <- S[orientate,]$CASI_Cg * (-1)
  
  S[orientate,]$Cbp_Co <- S[orientate,]$Cbp_Co * (-1)
  S[orientate,]$ASI_Co <- S[orientate,]$ASI_Co * (-1)
  S[orientate,]$EUR_Co <- S[orientate,]$EUR_Co * (-1)
  S[orientate,]$ME_Co <- S[orientate,]$ME_Co * (-1)
  S[orientate,]$CASI_Co <- S[orientate,]$CASI_Co * (-1)
  
  # Calculate delta S by comparing the absolute values of S_Cbp_Cg and S_Cbp_Co.
  # (ΔS < 0 means dominance of CG and ΔS > 0 means dominance of CO)
  S$delta_Cbp <- abs(S$Cbp_Cg) - abs(S$Cbp_Co)
  S$delta_ASI <- abs(S$ASI_Cg) - abs(S$ASI_Co)
  S$delta_EUR <- abs(S$EUR_Cg) - abs(S$EUR_Co)
  S$delta_ME <- abs(S$ME_Cg) - abs(S$ME_Co)
  S$delta_CASI <- abs(S$CASI_Cg) - abs(S$CASI_Co)
  head(S)
  
  # second loop
  population<-c("Cbp_Cg","ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg","Cbp_Co","ASI_Co","EUR_Co","ME_Co","CASI_Co") # set population
  for (j in population) {
    # count proportion of S < 0
    prop <- length(which(S[,which(colnames(S) == j)] < 0)) / length(S[,which(colnames(S) == j)])

    # binormal test for S<0, S>0
    b_test<-binom.test(c(length(which(S[,which(colnames(S) == j)] < 0)),length(which(S[,which(colnames(S) == j)] > 0))),p=0.5)
    if(b_test$p.value>=0.05){ # set the left label
      leftOut<-paste(round(prop,2)*100,"% ns",sep="")}else{
        if(b_test$p.value<0.001){leftOut<-paste(round(prop,2)*100,"% ***",sep="")}else{
          if(b_test$p.value<0.01){leftOut<-paste(round(prop,2)*100,"% **",sep="")}else{
            leftOut<-paste(round(prop,2),"*",sep=" ")}}
      }
    rightOut <- paste0(round(1-prop,2)*100,"%") # set the right label
    myMedian <- paste0("median=",round(median(S[,which(colnames(S) == j)]),2))
    
    S$gg <- S[,which(colnames(S) == j)]
    # plot distribution of S
    library(ggplot2)
    
    if(sub(".*_", "", j) == "Cg"){
      pop <- sub("_Cg","",j)
      p1 <- ggplot(S, aes(x=gg)) + 
        geom_histogram(aes(y=..density..),
                       color="black", fill="#FFF0F5", bins = 40) +
        xlim(c(-1, 1)) + ylim(c(0,2)) +
        geom_vline(aes(xintercept=0),
                   color="gray23", linetype="dashed", size=1) +
        geom_vline(aes(xintercept=median(gg)),
                   color="red", linetype="dashed", size=1) +
        labs(x =expression(S[Cg]), y = "Density", title = paste0(pop, " ", Tissue, " (n=", dim(S)[1], ")")) +
        annotate(geom="text", x=-1, y=1.9, label="Towards CG", hjust=0) +
        annotate(geom="text", x=1, y=1.9, label="Towards CO", hjust=1) +
        annotate(geom="text", x=-1, y=1.9*0.9, label=leftOut, hjust=0) +
        annotate(geom="text", x=1, y=1.9*0.9, label=rightOut, hjust=1) +
        annotate(geom="text", x=-1, y=1.9*0.8, label=myMedian, color="red", hjust=0) +
        theme_classic()
      assign(paste0("p_", i,"_", j,"_all"), p1)
    }else{
      pop <- sub("_Co","",j)
      p1 <- ggplot(S, aes(x=gg)) + 
        geom_histogram(aes(y=..density..),
                       color="black", fill="#F0FFFF", bins = 40) +
        xlim(c(-1, 1)) + ylim(c(0,2)) +
        geom_vline(aes(xintercept=0),
                   color="gray23", linetype="dashed", size=1) +
        geom_vline(aes(xintercept=median(gg)),
                   color="blue", linetype="dashed", size=1) +
        labs(x =expression(S[Co]), y = "Density", title = paste0(pop, " ", Tissue, " (n=", dim(S)[1], ")")) +
        annotate(geom="text", x=-1, y=1.9, label="Towards CG", hjust=0) +
        annotate(geom="text", x=1, y=1.9, label="Towards CO", hjust=1) +
        annotate(geom="text", x=-1, y=1.9*0.9, label=leftOut, hjust=0) +
        annotate(geom="text", x=1, y=1.9*0.9, label=rightOut, hjust=1) +
        annotate(geom="text", x=-1, y=1.9*0.8, label=myMedian, color="blue", hjust=0) +
        theme_classic()
      assign(paste0("p_", i,"_", j,"_all"), p1)
    }
  }
  
  # Third Loop
  # DE genes in different comparisons
  for (com in c("CGCO", "CGCR", "CRCO")){
    DE <- read.table(paste0("OutputData/DE_", com, "_", i, "_FDR.05.txt"), header = T, sep = "\t")
    
    # keep genes with CRCO FDR ≤ 0.05
    keep <- S$Gene %in% row.names(DE)
    S_DE <- S[keep,]
    
    # second loop for plot S distribution in DE genes
    population<-c("Cbp_Cg","ASI_Cg","EUR_Cg","ME_Cg","CASI_Cg","Cbp_Co","ASI_Co","EUR_Co","ME_Co","CASI_Co") # set population
    for (k in population) {
      # count proportion of S < 0
      prop <- length(which(S_DE[,which(colnames(S_DE) == k)] < 0)) / length(S_DE[,which(colnames(S_DE) == k)])
      
      # binormal test for S<0, S>0
      b_test<-binom.test(c(length(which(S_DE[,which(colnames(S_DE) == k)] < 0)),length(which(S_DE[,which(colnames(S_DE) == k)] > 0))),p=0.5)
      if(b_test$p.value>=0.05){ # set the left label
        leftOut<-paste(round(prop,2)*100,"% ns",sep="")}else{
          if(b_test$p.value<0.001){leftOut<-paste(round(prop,2)*100,"% ***",sep="")}else{
            if(b_test$p.value<0.01){leftOut<-paste(round(prop,2)*100,"% **",sep="")}else{
              leftOut<-paste(round(prop,2),"*",sep=" ")}}
        }
      rightOut <- paste0(round(1-prop,2)*100,"%") # set the right label
      myMedian <- paste0("median=",round(median(S_DE[,which(colnames(S_DE) == k)]),2))
      
      S_DE$gg <- S_DE[,which(colnames(S_DE) == k)]
      # plot distribution of S
      library(ggplot2)
      
      if(sub(".*_", "", k) == "Cg"){
        pop <- sub("_Cg","",k)
        p2 <- ggplot(S_DE, aes(x=gg)) + 
          geom_histogram(aes(y=..density..),
                         color="black", fill="#FFF0F5", bins = 40) +
          xlim(c(-1, 1)) + ylim(c(0,2)) +
          geom_vline(aes(xintercept=0),
                     color="gray23", linetype="dashed", size=1) +
          geom_vline(aes(xintercept=median(gg)),
                     color="red", linetype="dashed", size=1) +
          labs(x =expression(S[Cg]), y = "Density", title = paste0(pop," ",com, " DEs in ",i," (n=", dim(S_DE)[1], ")")) +
          annotate(geom="text", x=-1, y=1.9, label="Towards CG", hjust=0) +
          annotate(geom="text", x=1, y=1.9, label="Towards CO", hjust=1) +
          annotate(geom="text", x=-1, y=1.9*0.9, label=leftOut, hjust=0) +
          annotate(geom="text", x=1, y=1.9*0.9, label=rightOut, hjust=1) +
          annotate(geom="text", x=-1, y=1.9*0.8, label=myMedian, color="red", hjust=0) +
          theme_classic()
        assign(paste0("p_", i,"_", k,"_", com, "_DEs"), p2)
      }else{
        pop <- sub("_Co","",k)
        p2 <- ggplot(S_DE, aes(x=gg)) + 
          geom_histogram(aes(y=..density..),
                         color="black", fill="#F0FFFF", bins = 40) +
          xlim(c(-1, 1)) + ylim(c(0,2)) +
          geom_vline(aes(xintercept=0),
                     color="gray23", linetype="dashed", size=1) +
          geom_vline(aes(xintercept=median(gg)),
                     color="blue", linetype="dashed", size=1) +
          labs(x =expression(S[Co]), y = "Density", title = paste0(pop," ",com," DEs in ",i," (n=", dim(S_DE)[1], ")")) +
          annotate(geom="text", x=-1, y=1.9, label="Towards CG", hjust=0) +
          annotate(geom="text", x=1, y=1.9, label="Towards CO", hjust=1) +
          annotate(geom="text", x=-1, y=1.9*0.9, label=leftOut, hjust=0) +
          annotate(geom="text", x=1, y=1.9*0.9, label=rightOut, hjust=1) +
          annotate(geom="text", x=-1, y=1.9*0.8, label=myMedian, color="blue", hjust=0) +
          theme_classic()
        assign(paste0("p_", i,"_", k,"_",com,"_DEs"), p2)
      }
      
    }
    
    # Assign the S & S_DE datasets to new files with Tissue name.
    assign(paste0("S_",i), S)
    write.table(S, file = paste0("OutputData/CGCO_S_All_",Tissue,".txt"), quote = F, sep = "\t", row.names = F)
    
    assign(paste0("S_DE_",com,"_",i), S_DE)
    write.table(S_DE, file = paste0("OutputData/CGCO_S_",com,"_DEs_",Tissue,".txt"), quote = F, sep = "\t", row.names = F)
    
    ######################################################
    # Convergence index
    # for each transcript i, we used a convergence index Ci = ∆parents - ∆x / max(∆parents, ∆x)
    # ∆x stands for either ΔCo, ΔCg or Δsub
    # Ci thus ranges from -1 to 1, with positive values indicating more similar expression between the subgenomes
    # of C. bursa-pastoris than between parental species, and negative values indicating increased diffeRences
    # between subgenomes; the closer Ci to 0, the more similar are the expression patterns to parental species.
    
    S_DE$gg <- NULL # remove the colume specfic for ggplot
    library(MatrixGenerics)
    # Absolute differences between parents in expression.
    S_DE$C_delta_Parents <- abs(S_DE$mean_CG - S_DE$mean_CO)
    
    # Absolute differences between Subgenomes in expression.
    S_DE$C_delta_Cbp <- abs(S_DE$mean_Cbp_Cg - S_DE$mean_Cbp_Co) # Cbp all
    S_DE$C_delta_ASI <- abs(S_DE$mean_ASI_Cg - S_DE$mean_ASI_Co) # Cbp-ASI
    S_DE$C_delta_EUR <- abs(S_DE$mean_EUR_Cg - S_DE$mean_EUR_Co) # Cbp-EUR
    S_DE$C_delta_ME <- abs(S_DE$mean_ME_Cg - S_DE$mean_ME_Co) # Cbp-ME
    S_DE$C_delta_CASI <- abs(S_DE$mean_CASI_Cg - S_DE$mean_CASI_Co) # Cbp-CASI
    
    # each subgenome and the opposite parental species
    S_DE$C_delta_Cbp_cg <- abs(S_DE$mean_Cbp_Cg - S_DE$mean_CO) # Cbp all
    S_DE$C_delta_Cbp_co <- abs(S_DE$mean_Cbp_Co - S_DE$mean_CG)
    S_DE$C_delta_ASI_cg <- abs(S_DE$mean_ASI_Cg - S_DE$mean_CO) # ASI
    S_DE$C_delta_ASI_co <- abs(S_DE$mean_ASI_Co - S_DE$mean_CG)
    S_DE$C_delta_EUR_cg <- abs(S_DE$mean_EUR_Cg - S_DE$mean_CO) # EUR
    S_DE$C_delta_EUR_co <- abs(S_DE$mean_EUR_Co - S_DE$mean_CG)
    S_DE$C_delta_ME_cg <- abs(S_DE$mean_ME_Cg - S_DE$mean_CO) # ME
    S_DE$C_delta_ME_co <- abs(S_DE$mean_ME_Co - S_DE$mean_CG)
    S_DE$C_delta_CASI_cg <- abs(S_DE$mean_CASI_Cg - S_DE$mean_CO) # CASI
    S_DE$C_delta_CASI_co <- abs(S_DE$mean_CASI_Co - S_DE$mean_CG)
    
    # define columns with empty value.
    S_DE$C.Cbp <- NA # Cbp
    S_DE$C.Cbp_cg <- NA
    S_DE$C.Cbp_co <- NA
    S_DE$C.ASI <- NA # ASI
    S_DE$C.ASI_cg <- NA
    S_DE$C.ASI_co <- NA
    S_DE$C.EUR <- NA # EUR
    S_DE$C.EUR_cg <- NA
    S_DE$C.EUR_co <- NA
    S_DE$C.ME <- NA # ME
    S_DE$C.ME_cg <- NA
    S_DE$C.ME_co <- NA
    S_DE$C.CASI <- NA # CASI
    S_DE$C.CASI_cg <- NA
    S_DE$C.CASI_co <- NA
    
    for(a in c(1:dim(S_DE)[1])){
      S_DE[a,]$C.Cbp <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_Cbp) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_Cbp")]))
      S_DE[a,]$C.Cbp_cg <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_Cbp_cg) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_Cbp_cg")]))
      S_DE[a,]$C.Cbp_co <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_Cbp_co) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_Cbp_co")]))
      
      S_DE[a,]$C.ASI <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ASI) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ASI")]))
      S_DE[a,]$C.ASI_cg <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ASI_cg) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ASI_cg")]))
      S_DE[a,]$C.ASI_co <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ASI_co) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ASI_co")]))
      
      S_DE[a,]$C.EUR <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_EUR) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_EUR")]))
      S_DE[a,]$C.EUR_cg <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_EUR_cg) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_EUR_cg")]))
      S_DE[a,]$C.EUR_co <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_EUR_co) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_EUR_co")]))
      
      S_DE[a,]$C.ME <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ME) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ME")]))
      S_DE[a,]$C.ME_cg <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ME_cg) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ME_cg")]))
      S_DE[a,]$C.ME_co <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_ME_co) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_ME_co")]))
      
      S_DE[a,]$C.CASI <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_CASI) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_CASI")]))
      S_DE[a,]$C.CASI_cg <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_CASI_cg) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_CASI_cg")]))
      S_DE[a,]$C.CASI_co <- (S_DE[a,]$C_delta_Parents - S_DE[a,]$C_delta_CASI_co) / rowMaxs(as.matrix(S_DE[a, c("C_delta_Parents","C_delta_CASI_co")]))
    }
    
    # Assign the S & C DE datasets to new files with Tissue name.
    assign(paste0("SC_DE_",com,"_",i), S_DE)
    write.table(S_DE, file = paste0("OutputData/CGCO_SC_",com,"_DEs_",Tissue,".txt"), quote = F, sep = "\t", row.names = F)
    
  }
  
  
  
}

library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load

# Cbp in all populations
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,2)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_Cbp_Cg_all, vp = vplayout(1,1))
print(p_F_Cbp_Co_all, vp = vplayout(1,2))
print(p_L_Cbp_Cg_all, vp = vplayout(2,1))
print(p_L_Cbp_Co_all, vp = vplayout(2,2))
print(p_R_Cbp_Cg_all, vp = vplayout(3,1))
print(p_R_Cbp_Co_all, vp = vplayout(3,2))

# Cbp in each population, ASI, EUR, ME, CASI
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_ASI_Cg_all, vp = vplayout(1,1))
print(p_F_ASI_Co_all, vp = vplayout(1,2))
print(p_F_EUR_Cg_all, vp = vplayout(2,1))
print(p_F_EUR_Co_all, vp = vplayout(2,2))
print(p_F_ME_Cg_all, vp = vplayout(3,1))
print(p_F_ME_Co_all, vp = vplayout(3,2))
print(p_F_CASI_Cg_all, vp = vplayout(4,1))
print(p_F_CASI_Co_all, vp = vplayout(4,2))
print(p_L_ASI_Cg_all, vp = vplayout(1,3))
print(p_L_ASI_Co_all, vp = vplayout(1,4))
print(p_L_EUR_Cg_all, vp = vplayout(2,3))
print(p_L_EUR_Co_all, vp = vplayout(2,4))
print(p_L_ME_Cg_all, vp = vplayout(3,3))
print(p_L_ME_Co_all, vp = vplayout(3,4))
print(p_L_CASI_Cg_all, vp = vplayout(4,3))
print(p_L_CASI_Co_all, vp = vplayout(4,4))
print(p_R_ASI_Cg_all, vp = vplayout(1,5))
print(p_R_ASI_Co_all, vp = vplayout(1,6))
print(p_R_EUR_Cg_all, vp = vplayout(2,5))
print(p_R_EUR_Co_all, vp = vplayout(2,6))
print(p_R_ME_Cg_all, vp = vplayout(3,5))
print(p_R_ME_Co_all, vp = vplayout(3,6))
print(p_R_CASI_Cg_all, vp = vplayout(4,5))
print(p_R_CASI_Co_all, vp = vplayout(4,6))

# S_DE Cbp in all populations
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_Cbp_Cg_all, vp = vplayout(1,1))
print(p_F_Cbp_Co_all, vp = vplayout(1,2))
print(p_L_Cbp_Cg_all, vp = vplayout(1,3))
print(p_L_Cbp_Co_all, vp = vplayout(1,4))
print(p_R_Cbp_Cg_all, vp = vplayout(1,5))
print(p_R_Cbp_Co_all, vp = vplayout(1,6))
print(p_F_Cbp_Cg_CRCO_DEs, vp = vplayout(2,1))
print(p_F_Cbp_Co_CRCO_DEs, vp = vplayout(2,2))
print(p_L_Cbp_Cg_CRCO_DEs, vp = vplayout(2,3))
print(p_L_Cbp_Co_CRCO_DEs, vp = vplayout(2,4))
print(p_R_Cbp_Cg_CRCO_DEs, vp = vplayout(2,5))
print(p_R_Cbp_Co_CRCO_DEs, vp = vplayout(2,6))
print(p_F_Cbp_Cg_CGCR_DEs, vp = vplayout(3,1))
print(p_F_Cbp_Co_CGCR_DEs, vp = vplayout(3,2))
print(p_L_Cbp_Cg_CGCR_DEs, vp = vplayout(3,3))
print(p_L_Cbp_Co_CGCR_DEs, vp = vplayout(3,4))
print(p_R_Cbp_Cg_CGCR_DEs, vp = vplayout(3,5))
print(p_R_Cbp_Co_CGCR_DEs, vp = vplayout(3,6))
print(p_F_Cbp_Cg_CGCO_DEs, vp = vplayout(4,1))
print(p_F_Cbp_Co_CGCO_DEs, vp = vplayout(4,2))
print(p_L_Cbp_Cg_CGCO_DEs, vp = vplayout(4,3))
print(p_L_Cbp_Co_CGCO_DEs, vp = vplayout(4,4))
print(p_R_Cbp_Cg_CGCO_DEs, vp = vplayout(4,5))
print(p_R_Cbp_Co_CGCO_DEs, vp = vplayout(4,6))

#S_DE Cbp in each population, ASI, EUR, ME, CASI
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_ASI_Cg_CRCO_DEs, vp = vplayout(1,1))
print(p_F_ASI_Co_CRCO_DEs, vp = vplayout(1,2))
print(p_F_EUR_Cg_CRCO_DEs, vp = vplayout(2,1))
print(p_F_EUR_Co_CRCO_DEs, vp = vplayout(2,2))
print(p_F_ME_Cg_CRCO_DEs, vp = vplayout(3,1))
print(p_F_ME_Co_CRCO_DEs, vp = vplayout(3,2))
print(p_F_CASI_Cg_CRCO_DEs, vp = vplayout(4,1))
print(p_F_CASI_Co_CRCO_DEs, vp = vplayout(4,2))
print(p_L_ASI_Cg_CRCO_DEs, vp = vplayout(1,3))
print(p_L_ASI_Co_CRCO_DEs, vp = vplayout(1,4))
print(p_L_EUR_Cg_CRCO_DEs, vp = vplayout(2,3))
print(p_L_EUR_Co_CRCO_DEs, vp = vplayout(2,4))
print(p_L_ME_Cg_CRCO_DEs, vp = vplayout(3,3))
print(p_L_ME_Co_CRCO_DEs, vp = vplayout(3,4))
print(p_L_CASI_Cg_CRCO_DEs, vp = vplayout(4,3))
print(p_L_CASI_Co_CRCO_DEs, vp = vplayout(4,4))
print(p_R_ASI_Cg_CRCO_DEs, vp = vplayout(1,5))
print(p_R_ASI_Co_CRCO_DEs, vp = vplayout(1,6))
print(p_R_EUR_Cg_CRCO_DEs, vp = vplayout(2,5))
print(p_R_EUR_Co_CRCO_DEs, vp = vplayout(2,6))
print(p_R_ME_Cg_CRCO_DEs, vp = vplayout(3,5))
print(p_R_ME_Co_CRCO_DEs, vp = vplayout(3,6))
print(p_R_CASI_Cg_CRCO_DEs, vp = vplayout(4,5))
print(p_R_CASI_Co_CRCO_DEs, vp = vplayout(4,6))

#S_DE Cbp in each population, ASI, EUR, ME, CASI
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_ASI_Cg_CGCR_DEs, vp = vplayout(1,1))
print(p_F_ASI_Co_CGCR_DEs, vp = vplayout(1,2))
print(p_F_EUR_Cg_CGCR_DEs, vp = vplayout(2,1))
print(p_F_EUR_Co_CGCR_DEs, vp = vplayout(2,2))
print(p_F_ME_Cg_CGCR_DEs, vp = vplayout(3,1))
print(p_F_ME_Co_CGCR_DEs, vp = vplayout(3,2))
print(p_F_CASI_Cg_CGCR_DEs, vp = vplayout(4,1))
print(p_F_CASI_Co_CGCR_DEs, vp = vplayout(4,2))
print(p_L_ASI_Cg_CGCR_DEs, vp = vplayout(1,3))
print(p_L_ASI_Co_CGCR_DEs, vp = vplayout(1,4))
print(p_L_EUR_Cg_CGCR_DEs, vp = vplayout(2,3))
print(p_L_EUR_Co_CGCR_DEs, vp = vplayout(2,4))
print(p_L_ME_Cg_CGCR_DEs, vp = vplayout(3,3))
print(p_L_ME_Co_CGCR_DEs, vp = vplayout(3,4))
print(p_L_CASI_Cg_CGCR_DEs, vp = vplayout(4,3))
print(p_L_CASI_Co_CGCR_DEs, vp = vplayout(4,4))
print(p_R_ASI_Cg_CGCR_DEs, vp = vplayout(1,5))
print(p_R_ASI_Co_CGCR_DEs, vp = vplayout(1,6))
print(p_R_EUR_Cg_CGCR_DEs, vp = vplayout(2,5))
print(p_R_EUR_Co_CGCR_DEs, vp = vplayout(2,6))
print(p_R_ME_Cg_CGCR_DEs, vp = vplayout(3,5))
print(p_R_ME_Co_CGCR_DEs, vp = vplayout(3,6))
print(p_R_CASI_Cg_CGCR_DEs, vp = vplayout(4,5))
print(p_R_CASI_Co_CGCR_DEs, vp = vplayout(4,6))


#S_DE Cbp in each population, ASI, EUR, ME, CASI
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_F_ASI_Cg_CGCO_DEs, vp = vplayout(1,1))
print(p_F_ASI_Co_CGCO_DEs, vp = vplayout(1,2))
print(p_F_EUR_Cg_CGCO_DEs, vp = vplayout(2,1))
print(p_F_EUR_Co_CGCO_DEs, vp = vplayout(2,2))
print(p_F_ME_Cg_CGCO_DEs, vp = vplayout(3,1))
print(p_F_ME_Co_CGCO_DEs, vp = vplayout(3,2))
print(p_F_CASI_Cg_CGCO_DEs, vp = vplayout(4,1))
print(p_F_CASI_Co_CGCO_DEs, vp = vplayout(4,2))
print(p_L_ASI_Cg_CGCO_DEs, vp = vplayout(1,3))
print(p_L_ASI_Co_CGCO_DEs, vp = vplayout(1,4))
print(p_L_EUR_Cg_CGCO_DEs, vp = vplayout(2,3))
print(p_L_EUR_Co_CGCO_DEs, vp = vplayout(2,4))
print(p_L_ME_Cg_CGCO_DEs, vp = vplayout(3,3))
print(p_L_ME_Co_CGCO_DEs, vp = vplayout(3,4))
print(p_L_CASI_Cg_CGCO_DEs, vp = vplayout(4,3))
print(p_L_CASI_Co_CGCO_DEs, vp = vplayout(4,4))
print(p_R_ASI_Cg_CGCO_DEs, vp = vplayout(1,5))
print(p_R_ASI_Co_CGCO_DEs, vp = vplayout(1,6))
print(p_R_EUR_Cg_CGCO_DEs, vp = vplayout(2,5))
print(p_R_EUR_Co_CGCO_DEs, vp = vplayout(2,6))
print(p_R_ME_Cg_CGCO_DEs, vp = vplayout(3,5))
print(p_R_ME_Co_CGCO_DEs, vp = vplayout(3,6))
print(p_R_CASI_Cg_CGCO_DEs, vp = vplayout(4,5))
print(p_R_CASI_Co_CGCO_DEs, vp = vplayout(4,6))


################# Create a new dataframe of S for plot ###################
# S indices in different cases (All, phylogeny, MST, Phylogeny + MST)
CGCO_S_case <- rep(c(rep("All",3), rep("Phylogeny",3), rep("MST", 3), rep("Phylogeny + MST",3)),3) # 3 tissues
CGCO_S_index <- c(rep(c("delta_S","S_Co", "S_Cg"), 12)) # 4 cases, 3 tissues (F,L,R)
CGCO_S_tissue <- c(rep("Flower", 12), rep("Leaf", 12), rep("Root", 12))
CGCO_S_value <- c(median((S_F$Cbp_Co) + (S_F$Cbp_Cg)), # Flower all
             median(S_F$Cbp_Co),
             median(S_F$Cbp_Cg),
             median((S_DE_CRCO_F$Cbp_Co) + (S_DE_CRCO_F$Cbp_Cg)), # Flower CRCO DEs
             median(S_DE_CRCO_F$Cbp_Co),
             median(S_DE_CRCO_F$Cbp_Cg),
             median(((S_DE_CGCR_F$Cbp_Co) + (S_DE_CGCR_F$Cbp_Cg))), # Flower CGCR DEs
             median(S_DE_CGCR_F$Cbp_Co),
             median(S_DE_CGCR_F$Cbp_Cg),
             median(((S_DE_CGCO_F$Cbp_Co) + (S_DE_CGCO_F$Cbp_Cg))), # Flower CGCO DEs
             median(S_DE_CGCO_F$Cbp_Co),
             median(S_DE_CGCO_F$Cbp_Cg),
             median(((S_L$Cbp_Co) + (S_L$Cbp_Cg))), # Leaf all
             median(S_L$Cbp_Co),
             median(S_L$Cbp_Cg),
             median(((S_DE_CRCO_L$Cbp_Co) + (S_DE_CRCO_L$Cbp_Cg))), # Leaf CRCO DEs
             median(S_DE_CRCO_L$Cbp_Co),
             median(S_DE_CRCO_L$Cbp_Cg),
             median(((S_DE_CGCR_L$Cbp_Co) + (S_DE_CGCR_L$Cbp_Cg))), # Leaf CGCR DEs
             median(S_DE_CGCR_L$Cbp_Co),
             median(S_DE_CGCR_L$Cbp_Cg),
             median(((S_DE_CGCO_L$Cbp_Co) + (S_DE_CGCO_L$Cbp_Cg))), # Leaf CGCO DEs
             median(S_DE_CGCO_L$Cbp_Co),
             median(S_DE_CGCO_L$Cbp_Cg),
             median(((S_R$Cbp_Co) + (S_R$Cbp_Cg))), # Root all
             median(S_R$Cbp_Co),
             median(S_R$Cbp_Cg),
             median(((S_DE_CRCO_R$Cbp_Co) + (S_DE_CRCO_R$Cbp_Cg))), # Root CRCO DEs
             median(S_DE_CRCO_R$Cbp_Co),
             median(S_DE_CRCO_R$Cbp_Cg),
             median((median(S_DE_CGCR_R$Cbp_Co) + (S_DE_CGCR_R$Cbp_Cg))), # Root CGCR DEs
             median(S_DE_CGCR_R$Cbp_Co),
             median(S_DE_CGCR_R$Cbp_Cg),
             median((median(S_DE_CGCO_R$Cbp_Co) + (S_DE_CGCO_R$Cbp_Cg))), # Root CGCO DEs
             median(S_DE_CGCO_R$Cbp_Co),
             median(S_DE_CGCO_R$Cbp_Cg)
             )
# Combine all vectors to a dataframe
CGCO_S_plot <- cbind(CGCO_S_case, CGCO_S_index, CGCO_S_tissue, CGCO_S_value)
# change the format of the dataframe
CGCO_S_plot <- print.data.frame(data.frame(CGCO_S_plot), quote=FALSE)
str(CGCO_S_plot)
CGCO_S_plot$CGCO_S_index <- factor(CGCO_S_plot$CGCO_S_index, levels = c("S_Co","S_Cg","delta_S"))
CGCO_S_plot$CGCO_S_case <- factor(CGCO_S_plot$CGCO_S_case, levels = c("All","Phylogeny","MST","Phylogeny + MST"))
CGCO_S_plot$CGCO_S_tissue <- factor(CGCO_S_plot$CGCO_S_tissue, levels = c("Flower","Leaf","Root"))
CGCO_S_plot$CGCO_S_value <- as.numeric(CGCO_S_plot$CGCO_S_value)
str(CGCO_S_plot)
CGCO_S_plot$Comparison <- "CG2CO"
names(CGCO_S_plot) <- c("Case","Index","Tissue","Value","Comparison")

write.table(CGCO_S_plot, file = "OutputData/ggplot2PlotData for S indecies in CGCO.txt", col.names = T, row.names = F, sep = "\t")

# RUN S C indices of CRCO in loop.R
S_plot <- rbind(CGCO_S_plot,CRCO_S_plot,CGCR_S_plot)
write.table(S_plot, file = "OutputData/ggplot2PlotData for S indecies in AllKinds.txt", col.names = T, row.names = F, sep = "\t")

# Plot
library(ggplot2)
ggplot(CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  facet_grid( ~ Case) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

head(S_plot)
p1 <- ggplot(CGCO_S_plot[CGCO_S_plot$Case == "All",], aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.15, 0.15)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

p2 <- ggplot(CRCO_S_plot[CRCO_S_plot$Case == "All",], aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.15, 0.15)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

p3 <- ggplot(CGCR_S_plot[CGCR_S_plot$Case == "All",], aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.15, 0.15)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

grid.newpage()
library(ggpubr)
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right") # set the common legend



# S indices in different populations
S_index <- c(rep(c("delta_S","S_Co", "S_Cg"), 15)) # 5 Cbp subgroups (all + asi + eur + me + casi), 3 tissues (F,L,R)
S_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 3)
S_tissue <- c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15))
# Count median S value for each subpopulation. Median value
S_value <- c(abs(median(S_DE_F$Cbp_Co)) - abs(median(S_DE_F$Cbp_Cg)), # Flower all  ##### TMM with cgco 
                      median(S_DE_F$Cbp_Co),
                      median(S_DE_F$Cbp_Cg),
                      abs(median(S_DE_F$ASI_Co)) - abs(median(S_DE_F$ASI_Cg)), # Flower ASI
                      median(S_DE_F$ASI_Co),
                      median(S_DE_F$ASI_Cg),
                      abs(median(S_DE_F$EUR_Co)) - abs(median(S_DE_F$EUR_Cg)), # Flower EUR
                      median(S_DE_F$EUR_Co),
                      median(S_DE_F$EUR_Cg),
                      abs(median(S_DE_F$ME_Co)) - abs(median(S_DE_F$ME_Cg)), # Flower ME
                      median(S_DE_F$ME_Co),
                      median(S_DE_F$ME_Cg),
                      abs(median(S_DE_F$CASI_Co)) - abs(median(S_DE_F$CASI_Cg)), # Flower CASI
                      median(S_DE_F$CASI_Co),
                      median(S_DE_F$CASI_Cg),
                      abs(median(S_DE_L$Cbp_Co)) - abs(median(S_DE_L$Cbp_Cg)), # Leaf all  ##### Leaf
                      median(S_DE_L$Cbp_Co),
                      median(S_DE_L$Cbp_Cg),
                      abs(median(S_DE_L$ASI_Co)) - abs(median(S_DE_L$ASI_Cg)), # Leaf ASI
                      median(S_DE_L$ASI_Co),
                      median(S_DE_L$ASI_Cg),
                      abs(median(S_DE_L$EUR_Co)) - abs(median(S_DE_L$EUR_Cg)), # Leaf EUR
                      median(S_DE_L$EUR_Co),
                      median(S_DE_L$EUR_Cg),
                      abs(median(S_DE_L$ME_Co)) - abs(median(S_DE_L$ME_Cg)), # Leaf ME
                      median(S_DE_L$ME_Co),
                      median(S_DE_L$ME_Cg),
                      abs(median(S_DE_L$CASI_Co)) - abs(median(S_DE_L$CASI_Cg)), # Leaf CASI
                      median(S_DE_L$CASI_Co),
                      median(S_DE_L$CASI_Cg),
                      abs(median(S_DE_R$Cbp_Co)) - abs(median(S_DE_R$Cbp_Cg)), # Root all  ##### Root
                      median(S_DE_R$Cbp_Co),
                      median(S_DE_R$Cbp_Cg),
                      abs(median(S_DE_R$ASI_Co)) - abs(median(S_DE_R$ASI_Cg)), # Root ASI
                      median(S_DE_R$ASI_Co),
                      median(S_DE_R$ASI_Cg),
                      abs(median(S_DE_R$EUR_Co)) - abs(median(S_DE_R$EUR_Cg)), # Root EUR
                      median(S_DE_R$EUR_Co),
                      median(S_DE_R$EUR_Cg),
                      abs(median(S_DE_R$ME_Co)) - abs(median(S_DE_R$ME_Cg)), # Root ME
                      median(S_DE_R$ME_Co),
                      median(S_DE_R$ME_Cg),
                      abs(median(S_DE_R$CASI_Co)) - abs(median(S_DE_R$CASI_Cg)), # Root CASI
                      median(S_DE_R$CASI_Co),
                      median(S_DE_R$CASI_Cg)
)


# Combine all vectors to a dataframe
S_plot <- cbind(S_index, S_population, S_tissue, S_value)
# change the format of the dataframe
S_plot <- print.data.frame(data.frame(S_plot), quote=FALSE)
str(S_plot)
S_plot$S_index <- factor(S_plot$S_index, levels = c("S_Co","S_Cg","delta_S"))
S_plot$S_population <- factor(S_plot$S_population, levels = c("Cbp","Cbp_ASI","Cbp_EUR","Cbp_ME","Cbp_CASI"))
S_plot$S_tissue <- factor(S_plot$S_tissue, levels = c("Flower","Leaf","Root"))
S_plot$S_value <- as.numeric(S_plot$S_value)
str(S_plot)

# Plot
library(ggplot2)
p_s <- ggplot(S_plot, aes(x=S_tissue, y=S_value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=S_index, shape=S_index),size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.4, 0.4)) +
  scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
                     labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  facet_grid( ~ S_population) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()


###########   Convergence index   ###########
# Create a new dataframe of C for plot
C_index <- c(rep(c("Cbp","Cbp_co", "Cbp_cg"), 15)) # 5 Cbp subgroups (all + asi + eur + me + casi), 3 tissues (F,L,R)
C_population <- rep(c(rep("Cbp", 3), rep("Cbp_ASI", 3), rep("Cbp_EUR", 3), rep("Cbp_ME", 3), rep("Cbp_CASI", 3)), 3)
C_tissue <- c(rep("Flower", 15), rep("Leaf", 15), rep("Root", 15))
# Count proportion of C > 0
C_value <- c(nrow(SC_DE_F[SC_DE_F$C.Cbp > 0, ]) / nrow(SC_DE_F), # Flower all  ##### Flower
             nrow(SC_DE_F[SC_DE_F$C.Cbp_co > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.Cbp_cg > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.ASI > 0, ]) / nrow(SC_DE_F), # Flower ASI
             nrow(SC_DE_F[SC_DE_F$C.ASI_co > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.ASI_cg > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.EUR > 0, ]) / nrow(SC_DE_F), # Flower EUR
             nrow(SC_DE_F[SC_DE_F$C.EUR_co > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.EUR_cg > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.ME > 0, ]) / nrow(SC_DE_F), # Flower ME
             nrow(SC_DE_F[SC_DE_F$C.ME_co > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.ME_cg > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.CASI > 0, ]) / nrow(SC_DE_F), # Flower CASI
             nrow(SC_DE_F[SC_DE_F$C.CASI_co > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_F[SC_DE_F$C.CASI_cg > 0, ]) / nrow(SC_DE_F),
             nrow(SC_DE_L[SC_DE_L$C.Cbp > 0, ]) / nrow(SC_DE_L), # Leaf all  ##### Leaf
             nrow(SC_DE_L[SC_DE_L$C.Cbp_co > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.Cbp_cg > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.ASI > 0, ]) / nrow(SC_DE_L), # Leaf ASI
             nrow(SC_DE_L[SC_DE_L$C.ASI_co > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.ASI_cg > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.EUR > 0, ]) / nrow(SC_DE_L), # Leaf EUR
             nrow(SC_DE_L[SC_DE_L$C.EUR_co > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.EUR_cg > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.ME > 0, ]) / nrow(SC_DE_L), # Leaf ME
             nrow(SC_DE_L[SC_DE_L$C.ME_co > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.ME_cg > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.CASI > 0, ]) / nrow(SC_DE_L), # Leaf CASI
             nrow(SC_DE_L[SC_DE_L$C.CASI_co > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_L[SC_DE_L$C.CASI_cg > 0, ]) / nrow(SC_DE_L),
             nrow(SC_DE_R[SC_DE_R$C.Cbp > 0, ]) / nrow(SC_DE_R), # Root all  ##### Root
             nrow(SC_DE_R[SC_DE_R$C.Cbp_co > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.Cbp_cg > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.ASI > 0, ]) / nrow(SC_DE_R), # Root ASI
             nrow(SC_DE_R[SC_DE_R$C.ASI_co > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.ASI_cg > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.EUR > 0, ]) / nrow(SC_DE_R), # Root EUR
             nrow(SC_DE_R[SC_DE_R$C.EUR_co > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.EUR_cg > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.ME > 0, ]) / nrow(SC_DE_R), # Root ME
             nrow(SC_DE_R[SC_DE_R$C.ME_co > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.ME_cg > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.CASI > 0, ]) / nrow(SC_DE_R), # Root CASI
             nrow(SC_DE_R[SC_DE_R$C.CASI_co > 0, ]) / nrow(SC_DE_R),
             nrow(SC_DE_R[SC_DE_R$C.CASI_cg > 0, ]) / nrow(SC_DE_R)
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

p_c <- ggplot(C_plot, aes(x = C_tissue, y = C_value, colour= C_index)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey31") +
  geom_point(aes(shape=C_index), position=position_dodge(width=0.5), size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(3, 4, 8))+
  scale_color_manual(values=c("orange3","dodgerblue", "brown1"))+
  #scale_color_manual(values=c('#E69F00','#56B4E9',"tomato" ))+
  scale_y_continuous(limits=c(0.4, 0.8), breaks = seq(0.4, 0.8, by = 0.1)) +
  scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
                  labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  facet_grid(. ~ C_population) + 
  labs(x ="Tissues", y = expression("proportion of " * C[i] * " > 0") ) +
  theme_bw()

library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p_s, vp = vplayout(1,1))
print(p_c, vp = vplayout(2,1))






