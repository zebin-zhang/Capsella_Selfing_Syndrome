#  D index pairwise comparison among Cg, Cr, and Co.
## Set work directory and load files
setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
#### This file contains reads counts infromation of each gene in all samples
TMM <- read.table("InputData/Diploids_individual_TMM_FLR.txt", header=TRUE)
head(TMM)
dim(TMM)

### Define populations
CR_F <- c("CR1_F","CR2_F","CR3_F","CR4_F")
CG_F <- c("CG1_F","CG2_F","CG3_F","CG4_F")
CO_F <- c("CO1_F","CO2_F","CO3_F","CO4_F")

CR_L <- c("CR1_L","CR2_L","CR3_L","CR4_L")
CG_L <- c("CG1_L","CG2_L","CG3_L","CG4_L")
CO_L <- c("CO1_L","CO2_L","CO3_L","CO4_L")

CR_R <- c("CR1_R","CR2_R","CR3_R","CR4_R")
CG_R <- c("CG1_R","CG2_R","CG3_R","CG4_R")
CO_R <- c("CO1_R","CO2_R","CO3_R","CO4_R")

##############################
head(TMM)

TMM.Pop <- data.frame(Gene=row.names(TMM))
head(TMM.Pop)
# Use the MEAN value of each sample to represents as population value
# Flower of Cr, Cg, Co
TMM.Pop$CR_F <- rowMeans(TMM[,CR_F], na.rm = T)
TMM.Pop$CG_F <- rowMeans(TMM[,CG_F], na.rm = T)
TMM.Pop$CO_F <- rowMeans(TMM[,CO_F], na.rm = T)
# Leaf of Cr, Cg, Co
TMM.Pop$CR_L <- rowMeans(TMM[,CR_L], na.rm = T)
TMM.Pop$CG_L <- rowMeans(TMM[,CG_L], na.rm = T)
TMM.Pop$CO_L <- rowMeans(TMM[,CO_L], na.rm = T)
# Root of Cr, Cg, Co
TMM.Pop$CR_R <- rowMeans(TMM[,CR_R], na.rm = T)
TMM.Pop$CG_R <- rowMeans(TMM[,CG_R], na.rm = T)
TMM.Pop$CO_R <- rowMeans(TMM[,CO_R], na.rm = T)

row.names(TMM.Pop) <- as.character(TMM.Pop$Gene)
TMM.Pop$Gene <- NULL
head(TMM.Pop)
dim(TMM.Pop)
dim(na.omit(TMM.Pop))
str(TMM.Pop)

## Comparasion 
TMM.Pop.Comp <- TMM.Pop[,c("CG_F","CR_F","CO_F","CG_L","CR_L","CO_L","CG_R","CR_R","CO_R")]
head(TMM.Pop.Comp)
# First calcuate the absolute value between COCR and CGCR
### Flower
TMM.Pop.Comp$abs.F_CoCr <- abs(TMM.Pop.Comp[,"CO_F"] - TMM.Pop.Comp[,"CR_F"])
TMM.Pop.Comp$abs.F_CgCr <- abs(TMM.Pop.Comp[,"CG_F"] - TMM.Pop.Comp[,"CR_F"])
### Leaf
TMM.Pop.Comp$abs.L_CoCr <- abs(TMM.Pop.Comp[,"CO_L"] - TMM.Pop.Comp[,"CR_L"])
TMM.Pop.Comp$abs.L_CgCr <- abs(TMM.Pop.Comp[,"CG_L"] - TMM.Pop.Comp[,"CR_L"])
### Root
TMM.Pop.Comp$abs.R_CoCr <- abs(TMM.Pop.Comp[,"CO_R"] - TMM.Pop.Comp[,"CR_R"])
TMM.Pop.Comp$abs.R_CgCr <- abs(TMM.Pop.Comp[,"CG_R"] - TMM.Pop.Comp[,"CR_R"])

# calculate the discrepancy between cocr & cgcr
TMM.Pop.Comp$F.D.raw <- TMM.Pop.Comp$abs.F_CoCr - TMM.Pop.Comp$abs.F_CgCr
TMM.Pop.Comp$L.D.raw <- TMM.Pop.Comp$abs.L_CoCr - TMM.Pop.Comp$abs.L_CgCr
TMM.Pop.Comp$R.D.raw <- TMM.Pop.Comp$abs.R_CoCr - TMM.Pop.Comp$abs.R_CgCr

# install.packages("matrixStats")
library(matrixStats)
# Flower
TMM.Pop.Comp$F.max <- rowMaxs(as.matrix(TMM.Pop.Comp[,c("abs.F_CoCr","abs.F_CgCr")]), na.rm = TRUE)
# Leaf
TMM.Pop.Comp$L.max <- rowMaxs(as.matrix(TMM.Pop.Comp[,c("abs.L_CoCr","abs.L_CgCr")]), na.rm = TRUE)
# Root
TMM.Pop.Comp$R.max <- rowMaxs(as.matrix(TMM.Pop.Comp[,c("abs.R_CoCr","abs.R_CgCr")]), na.rm = TRUE)

# Perform the calculation
### Flower
TMM.Pop.Comp$F.D <- TMM.Pop.Comp$F.D.raw / TMM.Pop.Comp$F.max
TMM.Pop.Comp$L.D <- TMM.Pop.Comp$L.D.raw / TMM.Pop.Comp$L.max
TMM.Pop.Comp$R.D <- TMM.Pop.Comp$R.D.raw / TMM.Pop.Comp$R.max

head(TMM.Pop.Comp)
# Plot D index distribution in F L R.
# create new dataframe for F, L R
Value <- TMM.Pop.Comp$F.D  
D_F <- data.frame(Value)
D_F$Rep <- "Flower"
row.names(D_F) <- row.names(TMM.Pop.Comp)
head(D_F)
Value <- TMM.Pop.Comp$L.D
D_L <- data.frame(Value)
D_L$Rep <- "Leaf"
row.names(D_L) <- row.names(TMM.Pop.Comp)

Value <- TMM.Pop.Comp$R.D
D_R <- data.frame(Value)
D_R$Rep <- "Root"
row.names(D_R) <- row.names(TMM.Pop.Comp)

#####################################################################################################
# Analysis on UPPMAX
# 
# listOfDataFrames <- vector(mode = "list", length = 10000)  # Create a new list for all simulated data
# 
# for (i in 1:10000) {                    # set 10,000 simulations
#   set.seed(i)                           # set seed for each simulation for obtain certain values 
#   
#   TMM_random <- as.matrix(TMM[,1:12])
#   
#   for(j in 1:dim(TMM_random)[1]) {
#     TMM_random[j,] <- sample(as.matrix(TMM), 12)
#   }
#   
#   TMM_random.Pop <- data.frame(Gene=row.names(TMM_random))
#   # Use the MEAN value of each sample to represents as population value
#   # Flower of Cr, Cg, Co
#   TMM_random.Pop$CR_F <- rowMeans(TMM_random[,c(1:4)], na.rm = T)
#   TMM_random.Pop$CG_F <- rowMeans(TMM_random[,c(5:8)], na.rm = T)
#   TMM_random.Pop$CO_F <- rowMeans(TMM_random[,c(9:12)], na.rm = T)
# 
#   row.names(TMM_random.Pop) <- as.character(TMM_random.Pop$Gene)
#   TMM_random.Pop$Gene <- NULL
# 
#   ## Comparasion 
#   TMM_random.Pop.Comp <- TMM_random.Pop[,c("CG_F","CR_F","CO_F")]
#   head(TMM_random.Pop.Comp)
#   # First calcuate the absolute value between COCR and CGCR
#   ### Flower
#   TMM_random.Pop.Comp$abs.CoCr <- abs(TMM_random.Pop.Comp[,"CO_F"] - TMM_random.Pop.Comp[,"CR_F"])
#   TMM_random.Pop.Comp$abs.CgCr <- abs(TMM_random.Pop.Comp[,"CG_F"] - TMM_random.Pop.Comp[,"CR_F"])
#   # calculate the discrepancy between cocr & cgcr
#   TMM_random.Pop.Comp$D.raw <- TMM_random.Pop.Comp$abs.CoCr - TMM_random.Pop.Comp$abs.CgCr
# 
#   # install.packages("matrixStats")
#   library(matrixStats)
#   # Flower
#   TMM_random.Pop.Comp$max <- rowMaxs(as.matrix(TMM_random.Pop.Comp[,c("abs.CoCr","abs.CgCr")]), na.rm = TRUE)
#   # Perform the calculation
#   ### Flower
#   TMM_random.Pop.Comp$D <- TMM_random.Pop.Comp$D.raw / TMM_random.Pop.Comp$max
# 
#   if (i == 1){
#     Value <- TMM_random.Pop.Comp$D
#     Random <- data.frame(Value)
#     Random$Rep <- paste0("rep",i)
#     dim(Random)
#     listOfDataFrames[[i]] <- data.frame(Random)
#     Median <- median(Value)            # Calculate the median value for each simulation
#   }else{
#     Value <- TMM_random.Pop.Comp$D
#     App <- data.frame(Value)
#     App$Rep <- paste0("rep",i)
#     listOfDataFrames[[i]] <- data.frame(App)
#     Median[i] <- median(Value)
#   }
# }
# 
# Random <- do.call("rbind", listOfDataFrames) # combine all simulation by do.call
# 
# Random$Rep <- factor(Random$Rep)      # as factor
# Random_sd <- sd(Random$Value)
# greyScale <- rep(gray.colors(10, start = 0, end = 0.7), 1000)     # scale grey color to 10 scales from dark to light, repeat 100 times.
# 
# write.table(Random, file = "OutputData/Random_Simulation_D.txt", quote = F, sep = "\t", row.names = F, col.names = T)
# write.table(Median, file = "OutputData/Random_Simulation_D_Median.txt", quote = F, sep = "\t")
# getwd()
# Random <- read.table("OutputData/Random_Simulation_D.txt", header = T)
# head(Random)
# 
# density(Random[Random$Value > -0.7 & Random$Value < -0.3,]$Value)
# 
# 
# library(ggplot2)
# greyScale <- rep(gray.colors(10, start = 0, end = 0.7), 100)     # scale grey color to 10 scales from dark to light, repeat 100 times.
# 
# # Plot simulation result
# p <- ggplot(Random, aes(x=Value, color = Rep))+
#   geom_density() + 
#   scale_color_manual(values=greyScale) +   
#   xlim(-1,1) + ylim(0,1) +
#   theme_classic() +                                 # set theme to classic
#   theme(legend.position = "none") +                   # remove legend
#   labs( #title = "Fuel economy declines as weight increases",
#     #subtitle = "(1973-74)",
#     #caption = "Data from the 1974 Motor Trend US magazine.",
#     #tag = "D",
#     x = " ",
#     y = "Density"
#     #colour = "Gears"
#   )
# 
# d <- ggplot_build(p)
# 
# head(d$data) # check the density value
# Extract <- as.data.frame(d$data) # extract value
# head(Extract) 
# # Extract$density: density value of y axis
# # Extract$group: density groups, from 1 to 1000, 1000random simulations
# 
# # Extract values of each simulation within range of -0.7 : -0.3
# MyNegative <- Extract[Extract$x > -0.7 & Extract$x < -0.3,]
# head(MyNegative)
# calculate the minimum value of each simulation
# for(i in 1:1000){
#   if (i == 1){
#     MyMin <- min(MyNegative[MyNegative$group == i,]$density) # return the minimum density value
#     NegativeMin <- MyNegative[((MyNegative$group == i) & (round(MyNegative$density, 7) == round(MyMin,7))),]$x #return the x value, round is important in here
#   }else{
#     MyMin[i] <- min(MyNegative[MyNegative$group == i,]$density)
#     NegativeMin[i] <- MyNegative[((MyNegative$group == i) & (round(MyNegative$density, 7) == round(MyMin[i],7))),]$x #return the x value
#   }
# }
# NegativeMedian <- median(NegativeMin)
# 
# # Extract values of each simulation within range of 0.3 : 0.7
# MyPositive <- Extract[Extract$x > 0.3 & Extract$x < 0.7,]
# head(MyPositive)
# # calculate the minimum value of each simulation
# for(i in 1:1000){
#   if (i == 1){
#     MyMin <- min(MyPositive[MyPositive$group == i,]$density) # return the minimum density value
#     PositiveMin <- MyPositive[((MyPositive$group == i) & (round(MyPositive$density, 7) == round(MyMin,7))),]$x #return the x value, round is important in here
#   }else{
#     MyMin[i] <- min(MyPositive[MyPositive$group == i,]$density)
#     PositiveMin[i] <- MyPositive[((MyPositive$group == i) & (round(MyPositive$density, 7) == round(MyMin[i],7))),]$x #return the x value
#   }
# }
# PositiveMedian <- median(PositiveMin)
# 
# # Calculate percentage of positive, near zero, and negative
# head(head(Random)) 
# perPos <- dim(Random[Random$Value > PositiveMedian,])[1] /  dim(Random)[1]
# perZero <- dim(Random[(Random$Value <= PositiveMedian) & (Random$Value >= NegativeMedian),])[1] / dim(Random)[1]
# perNeg <- dim(Random[Random$Value < NegativeMedian,])[1] /  dim(Random)[1]
# 
# # plot simulation result
# 
# pdf("~/Desktop/D_10000_Simulations.pdf", width=5, height=4)
# pSimulation <- ggplot(Random, aes(x=Value, color = Rep))+
#   geom_density() + 
#   scale_color_manual(values=greyScale) +   
#   xlim(-1,1) + ylim(0,1) +
#   annotate(geom="text", x=-1, y=1, label=expression(D[italic(COCR)]~"<"~D[italic(CGCR)]), hjust=0) +
#   annotate(geom="text", x=-1, y=1*0.9, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
#   annotate(geom="text", x=0, y=1, label=expression(D[italic(COCR)]~"="~D[italic(CGCR)])) +
#   annotate(geom="text", x=0, y=1*0.9, label=paste(round(perZero, 4) * 100, "%")) +
#   annotate(geom="text", x=1, y=1, label=expression(D[italic(COCR)]~">"~D[italic(CGCR)]), hjust=1) +
#   annotate(geom="text", x=1, y=1*0.9, label=paste(round(perPos, 4) * 100, "%"), hjust=1) +
#   geom_vline(aes(xintercept=PositiveMedian), color="gray34", linetype="dashed", size=1) +
#   geom_vline(aes(xintercept=NegativeMedian), color="gray34", linetype="dashed", size=1) +
#   theme_classic() +                                 # set theme to classic
#   theme(legend.position = "none") +                   # remove legend
#   labs( #title = "Fuel economy declines as weight increases",
#     #subtitle = "(1973-74)",
#     #caption = "Data from the 1974 Motor Trend US magazine.",
#     #tag = "A",
#     x = " ",
#     y = "Density"
#     #colour = "Gears"
#   )
# 
# #dev.off()

library(ggplot2)
library(cowplot)

# set function for notin
'%!in%' <- function(x,y)!('%in%'(x,y))
################################### 1
tissue<-c("F","L","R") # set tissue
# Loop start
for(i in tissue){
  print(i)
  # set a sub loop for viable name
  if (i == "F"){
    D <- D_F
  }else{
    if (i == "L"){
      D <- D_L
    }else{
      D <- D_R
    }
  }
  # For nonDE
  CGCO <- read.csv(paste("OutputData/DE_", i, "_CGCO.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CGCO <- CGCO[CGCO$DE == "Not DE",]
  CGCR <- read.csv(paste("OutputData/DE_", i, "_CGCR.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CGCR <- CGCR[CGCR$DE == "Not DE",]
  CRCO <- read.csv(paste("OutputData/DE_", i, "_CRCO.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CRCO <- CRCO[CRCO$DE == "Not DE",]
  nonKeep <- (row.names(CGCO) %in% row.names(CGCR)) & (row.names(CGCO) %in% row.names(CRCO))
  nonDE_5 <- CGCO[nonKeep,] 
  CGCO <- read.csv(paste("OutputData/DE_", i, "_CGCO_FDR.01.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CGCO <- CGCO[CGCO$DE == "Not DE",]
  CGCR <- read.csv(paste("OutputData/DE_", i, "_CGCR_FDR.01.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CGCR <- CGCR[CGCR$DE == "Not DE",]
  CRCO <- read.csv(paste("OutputData/DE_", i, "_CRCO_FDR.01.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  CRCO <- CRCO[CRCO$DE == "Not DE",]
  nonKeep <- (row.names(CGCO) %in% row.names(CGCR)) & (row.names(CGCO) %in% row.names(CRCO))
  nonDE_1 <- CGCO[nonKeep,] 

  Density <- with(density(D$Value), data.frame(x, y))
  ########################################## 2
  # Seconde loop for All, DE-CGCR/CGCO,CRCO
  for (com in c("All", "FDR5", "FDR1", "CGCR", "CRCO", "CGCO")){
    Mylim <- 1.6
    if (com == "All"){
      Density <- with(density(D$Value, from = -1, to = 1), data.frame(x, y))
      N <- dim(D)[1]
      Pos <- round(median(D[D$Value > 0,]$Value), 2)
      Neg <- round(median(D[D$Value < 0,]$Value), 2)
      perPos <- dim(D[D$Value > 0,])[1] /  dim(D)[1]
      #perZero <- dim(D[(D$Value <= PositiveMedian) & (D$Value >= NegativeMedian),])[1] /  dim(D)[1]
      perNeg <- dim(D[D$Value < 0,])[1] /  dim(D)[1]
      de <- "all genes"
      MyMedian <- round(median(D$Value), 2)
    }else{
        if (com == "FDR5"){
          keep <- row.names(D) %!in% row.names(nonDE_5) # not in function
          subD <- D[keep,] 
          Density <- with(density(subD$Value, from = -1, to = 1), data.frame(x, y))
          N <- dim(subD)[1]
          Pos <- round(median(subD[subD$Value > 0,]$Value), 2)
          Neg <- round(median(subD[subD$Value < 0,]$Value), 2)
          perPos <- dim(subD[subD$Value > 0,])[1] /  dim(subD)[1]
          #perZero <- dim(subD[(subD$Value <= PositiveMedian) & (subD$Value >= NegativeMedian),])[1] /  dim(subD)[1]
          perNeg <- dim(subD[subD$Value < 0,])[1] /  dim(subD)[1]
          de <- "DE genes FDR < 0.05"
          MyMedian <- round(median(subD$Value), 2)
        }else{
          if (com == "FDR1"){
            keep <- row.names(D) %!in% row.names(nonDE_1) # not in function
            subD <- D[keep,] 
            Density <- with(density(subD$Value, from = -1, to = 1), data.frame(x, y))
            N <- dim(subD)[1]
            Pos <- round(median(subD[subD$Value > 0,]$Value), 2)
            Neg <- round(median(subD[subD$Value < 0,]$Value), 2)
            perPos <- dim(subD[subD$Value > 0,])[1] /  dim(subD)[1]
            #perZero <- dim(subD[(subD$Value <= PositiveMedian) & (subD$Value >= NegativeMedian),])[1] /  dim(subD)[1]
            perNeg <- dim(subD[subD$Value < 0,])[1] /  dim(subD)[1]
            de <- "DE genes FDR < 0.01"
            MyMedian <- round(median(subD$Value), 2)
          }else{
            DE <- read.csv(paste("OutputData/DE_", i, "_", com, "_FDR.05.txt", sep = ""), header = T, row.names = 1, sep = "\t")
            keep <- row.names(D) %in% row.names(DE)
            subD <- D[keep,] 
            MyN <- dim(subD)[1]
            assign(paste0("N_",com,"_",i), MyN)
            Density <- with(density(subD$Value, from = -1, to = 1), data.frame(x, y))
            N <- dim(subD)[1]
            Pos <- round(median(subD[subD$Value > 0,]$Value), 2)
            Neg <- round(median(subD[subD$Value < 0,]$Value), 2)
            perPos <- dim(subD[subD$Value > 0,])[1] /  dim(subD)[1]
            perNeg <- dim(subD[subD$Value < 0,])[1] /  dim(subD)[1]
            de <- "DE genes"
            MyMedian <- round(median(subD$Value), 2)
            write.table(subD, file = paste("OutputData/DE_", i, "_", com, "_FDR.05_D_index.txt", sep = ""), quote = F, sep = "\t")
          }
        }
    } # end if else
    assign(paste0(i, "_", com, "_Median"), MyMedian)
    assign(paste0(i, "_", com, "_Pos"), Pos)
    assign(paste0(i, "_", com, "_Neg"), Neg)
    assign(paste0(i, "_", com, "_perPos"), perPos)
    assign(paste0(i, "_", com, "_perNeg"), perNeg)
    assign(paste0(i, "_", com, "_Density"), Density)
    # Plot out
    if (i == "F"){
      p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
        geom_line(size=1 )+
        geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = "#4dbfea", size=1) + # light blue
        #geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
        geom_area(mapping = aes(x = ifelse(x < 0 , x, 0)), fill = "gray70") + 
        #geom_vline(aes(xintercept = Pos), color="gray34", linetype="dashed", size=1) +
        #geom_vline(aes(xintercept = 0), color="gray", linetype="dashed", size=1) +
        #geom_vline(aes(xintercept = Neg), color="red1", linetype="dashed", size=1) + 
        xlim(-1,1) + 
        ylim(0,Mylim) +
        annotate(geom="text", x=0, y=Mylim*0.9, label=paste("overall median =",MyMedian)) +
        annotate(geom="text", x=-1, y=Mylim, label= "D < 0", hjust=0) +
        annotate(geom="text", x=-1, y=Mylim*0.9, label=paste("median =",Neg), hjust=0) +
        annotate(geom="text", x=-1, y=Mylim*0.8, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
        annotate(geom="text", x=1, y=Mylim, label= "D > 0", hjust=1) +
        annotate(geom="text", x=1, y=Mylim*0.9, label=paste("median =",Pos), hjust=1) + 
        annotate(geom="text", x=1, y=Mylim*0.8, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
        theme_classic() +
        labs(title = paste("Flowers in", de,"(n=", N, ")"),
             x = "D index",
             y = "Density"
        ) 
      ggdraw(p)
      assign(paste0(i, "_", com, "_p"), p)
    }else{
      if (i == "L"){
        p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
          geom_line(size=1 )+
          geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = "#fae57d", size=1) + 
          #geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
          geom_area(mapping = aes(x = ifelse(x < 0 , x, 0)), fill = "gray90") + #light green
          #geom_vline(aes(xintercept = Pos), color="gray34", linetype="dashed", size=1) +
          #geom_vline(aes(xintercept = 0), color="gray", linetype="dashed", size=1) +
          #geom_vline(aes(xintercept = Neg), color="forestgreen", linetype="dashed", size=1) + 
          xlim(-1,1) + 
          ylim(0,Mylim) +
          annotate(geom="text", x=0, y=Mylim*0.9, label=paste("overall median =",MyMedian)) +
          annotate(geom="text", x=-1, y=Mylim, label= "D < 0", hjust=0) +
          annotate(geom="text", x=-1, y=Mylim*0.9, label=paste("median =",Neg), hjust=0) +
          annotate(geom="text", x=-1, y=Mylim*0.8, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
          annotate(geom="text", x=1, y=Mylim, label= "D > 0", hjust=1) +
          annotate(geom="text", x=1, y=Mylim*0.9, label=paste("median =",Pos), hjust=1) + 
          annotate(geom="text", x=1, y=Mylim*0.8, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
          theme_classic() +
          labs(title = paste("Leaves in", de,"(n=", N, ")"),
               x = "D index",
               y = "Density"
          ) 
        ggdraw(p)  
        assign(paste0(i, "_", com, "_p"), p)
      }else{
        p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
          geom_line(size=1 )+
          geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = "#fae57d", size=1) + 
          #geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
          geom_area(mapping = aes(x = ifelse(x < 0 , x, 0)), fill = "gray90") + #light brown
          #geom_vline(aes(xintercept = Pos), color="gray34", linetype="dashed", size=1) +
          #geom_vline(aes(xintercept = 0), color="gray", linetype="dashed", size=1) +
          #geom_vline(aes(xintercept = Neg), color="brown", linetype="dashed", size=1) + 
          xlim(-1,1) + 
          ylim(0,Mylim) +
          annotate(geom="text", x=0, y=Mylim*0.9, label=paste("overall median =",MyMedian)) +
          annotate(geom="text", x=-1, y=Mylim, label= "D < 0", hjust=0) +
          annotate(geom="text", x=-1, y=Mylim*0.9, label=paste("median =",Neg), hjust=0) +
          annotate(geom="text", x=-1, y=Mylim*0.8, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
          annotate(geom="text", x=1, y=Mylim, label= "D > 0", hjust=1) +
          annotate(geom="text", x=1, y=Mylim*0.9, label=paste("median =",Pos), hjust=1) + 
          annotate(geom="text", x=1, y=Mylim*0.8, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
          theme_classic() +
          labs(title = paste("Roots in", de,"(n=", N, ")"),
               x = "D index",
               y = "Density"
          ) 
        ggdraw(p)
        assign(paste0(i, "_", com, "_p"), p)
      }
    } # end if else
  } # end loop 2
} # end loop 1

library(ggpubr)
library("grid")
library(ggthemes) # Load
ggarrange(F_FDR1_p +
            geom_vline(aes(xintercept = F_FDR1_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = F_FDR1_Pos), color="#0067b2", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = F_FDR1_Neg), color="gray34", linetype="dashed", size=1),
          L_FDR1_p +
            geom_vline(aes(xintercept = L_FDR1_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = L_FDR1_Pos), color="#e8a628", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = L_FDR1_Neg), color="gray34", linetype="dashed", size=1),
          R_FDR1_p +
            geom_vline(aes(xintercept = R_FDR1_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = R_FDR1_Pos), color="#e8a628", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = R_FDR1_Neg), color="gray34", linetype="dashed", size=1),
          F_FDR5_p+
            geom_vline(aes(xintercept = F_FDR5_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = F_FDR5_Pos), color="#0067b2", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = F_FDR5_Neg), color="gray34", linetype="dashed", size=1),
          L_FDR5_p+
            geom_vline(aes(xintercept = L_FDR5_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = L_FDR5_Pos), color="#e8a628", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = L_FDR5_Neg), color="gray34", linetype="dashed", size=1),
          R_FDR5_p+
            geom_vline(aes(xintercept = R_FDR5_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = R_FDR5_Pos), color="#e8a628", linetype="dashed", size=1) + 
            geom_vline(aes(xintercept = R_FDR5_Neg), color="gray34", linetype="dashed", size=1),
          F_All_p + 
            geom_vline(aes(xintercept = F_All_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = F_All_Pos), color="#0067b2", linetype="dashed", size=1) + # vline can't add in loop
            geom_vline(aes(xintercept = F_All_Neg), color="gray34", linetype="dashed", size=1),
          L_All_p + 
            geom_vline(aes(xintercept = L_All_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = L_All_Pos), color="#e8a628", linetype="dashed", size=1) + # vline can't add in loop
            geom_vline(aes(xintercept = L_All_Neg), color="gray34", linetype="dashed", size=1),
          R_All_p + 
            geom_vline(aes(xintercept = R_All_Median), color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept = R_All_Pos), color="#e8a628", linetype="dashed", size=1) + # vline can't add in loop
            geom_vline(aes(xintercept = R_All_Neg), color="gray34", linetype="dashed", size=1),
          #labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3, ncol = 3)


pdf("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/Capsella_PDF/New_D_index_10000_simulations_20x16.pdf", width=20, height=16)
# Cbp in all populations
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,4)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(pSimulation, vp = vplayout(1,1))
print(F_All_p, vp = vplayout(2,1))
print(F_CGCR_p, vp = vplayout(2,2))
print(F_CRCO_p, vp = vplayout(2,3))
print(F_CGCO_p, vp = vplayout(2,4))
print(L_All_p, vp = vplayout(3,1))
print(L_CGCR_p, vp = vplayout(3,2))
print(L_CRCO_p, vp = vplayout(3,3))
print(L_CGCO_p, vp = vplayout(3,4))
print(R_All_p, vp = vplayout(4,1))
print(R_CGCR_p, vp = vplayout(4,2))
print(R_CRCO_p, vp = vplayout(4,3))
print(R_CGCO_p, vp = vplayout(4,4))
dev.off()


head(F_All_Density)

# All
D.Density <- cbind(F_All_Density, L_All_Density$y, R_All_Density$y)
names(D.Density) <- c("D","Flowers","Leaves","Roots")
head(D.Density)


Mylim <- 1
Flim <- 3
Llim <- 2
Rlim <- 1
p_all <- ggplot(data = D.Density) +
  geom_line(mapping = aes(x = D, y = Flowers +2),  size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, NA), y = Flowers + 2), fill = "#e55b7e") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Flowers + 2), fill = "#fee2e3") + 
  geom_area(mapping = aes(x = D, y = 2), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Leaves + 1), size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Leaves + 1), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Leaves + 1), fill = "#e6d6c5") + 
  geom_area(mapping = aes(x = D, y = 1), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Roots), size=1 ) +
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Roots), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Roots), fill = "#e6d6c5") + 
  geom_segment(aes(x = F_All_Median, y = 2, xend = F_All_Median, yend = 3), color = "gray", linetype = "dashed", size=1) + # Flowers
  geom_segment(aes(x = F_All_Pos, y = 2, xend = F_All_Pos, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = F_All_Neg, y = 2, xend = F_All_Neg, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0.5, y=Flim *1.1, label= "D > 0") +
  annotate(geom="text", x=-0.5, y=Flim *1.1, label= "D < 0") +
  annotate(geom="text", x=0, y=Flim, label=paste("overall median =", F_All_Median)) +
  annotate(geom="text", x=-1, y=Flim, label=paste("median =", F_All_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Flim -0.2, label=paste(round(F_All_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Flim, label=paste("median =",F_All_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Flim -0.2, label=paste(round(F_All_perPos, 4) * 100, "%"), hjust=1) +
  geom_segment(aes(x = L_All_Median, y = 1, xend = L_All_Median, yend = 2), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = L_All_Pos, y = 1, xend = L_All_Pos, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = L_All_Neg, y = 1, xend = L_All_Neg, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Llim -0.1, label=paste("overall median =", L_All_Median)) +
  annotate(geom="text", x=-1, y=Llim -0.1, label=paste("median =", L_All_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Llim -0.3, label=paste(round(L_All_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Llim -0.1, label=paste("median =",L_All_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Llim -0.3, label=paste(round(L_All_perPos, 4) * 100, "%"), hjust=1) +   
  geom_segment(aes(x = R_All_Median, y = 0, xend = R_All_Median, yend = 1), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = R_All_Pos, y = 0, xend = R_All_Pos, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = R_All_Neg, y = 0, xend = R_All_Neg, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Rlim -0.1, label=paste("overall median =", R_All_Median)) +
  annotate(geom="text", x=-1, y=Rlim -0.1, label=paste("median =", R_All_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Rlim -0.3, label=paste(round(R_All_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Rlim -0.1, label=paste("median =",R_All_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Rlim -0.3, label=paste(round(R_All_perPos, 4) * 100, "%"), hjust=1) +   
  xlim(-1,1) + 
  #ylim(0,3.4) +
  scale_y_continuous(limits = c(0,3.4), breaks=c(0,1,2,3), labels = c("","Roots","Leaves","Flowers")) +
  theme_classic() +
  labs(title = "All genes",
       x = "D index",
       y = "Density"
  ) 
  
FDR1Density <- cbind(F_FDR1_Density, L_FDR1_Density$y, R_FDR1_Density$y)
names(FDR1Density) <- c("D","Flowers","Leaves","Roots")
head(FDR1Density)

p_fdr1 <- ggplot(data = FDR1Density) +
  geom_line(mapping = aes(x = D, y = Flowers +2),  size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, NA), y = Flowers + 2), fill = "#e55b7e") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Flowers + 2), fill = "#fee2e3") + 
  geom_area(mapping = aes(x = D, y = 2), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Leaves + 1), size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Leaves + 1), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Leaves + 1), fill = "#e6d6c5") + 
  geom_area(mapping = aes(x = D, y = 1), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Roots), size=1 ) +
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Roots), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Roots), fill = "#e6d6c5") + 
  geom_segment(aes(x = F_FDR1_Median, y = 2, xend = F_FDR1_Median, yend = 3), color = "gray", linetype = "dashed", size=1) + # Flowers
  geom_segment(aes(x = F_FDR1_Pos, y = 2, xend = F_FDR1_Pos, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = F_FDR1_Neg, y = 2, xend = F_FDR1_Neg, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0.5, y=Flim *1.1, label= "D > 0") +
  annotate(geom="text", x=-0.5, y=Flim *1.1, label= "D < 0") +
  annotate(geom="text", x=0, y=Flim, label=paste("overall median =", F_FDR1_Median)) +
  annotate(geom="text", x=-1, y=Flim, label=paste("median =", F_FDR1_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Flim -0.2, label=paste(round(F_FDR1_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Flim, label=paste("median =",F_FDR1_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Flim -0.2, label=paste(round(F_FDR1_perPos, 4) * 100, "%"), hjust=1) +
  geom_segment(aes(x = L_FDR1_Median, y = 1, xend = L_FDR1_Median, yend = 2), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = L_FDR1_Pos, y = 1, xend = L_FDR1_Pos, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = L_FDR1_Neg, y = 1, xend = L_FDR1_Neg, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Llim -0.1, label=paste("overall median =", L_FDR1_Median)) +
  annotate(geom="text", x=-1, y=Llim -0.1, label=paste("median =", L_FDR1_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Llim -0.3, label=paste(round(L_FDR1_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Llim -0.1, label=paste("median =",L_FDR1_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Llim -0.3, label=paste(round(L_FDR1_perPos, 4) * 100, "%"), hjust=1) +   
  geom_segment(aes(x = R_FDR1_Median, y = 0, xend = R_FDR1_Median, yend = 1), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = R_FDR1_Pos, y = 0, xend = R_FDR1_Pos, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = R_FDR1_Neg, y = 0, xend = R_FDR1_Neg, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Rlim -0.1, label=paste("overall median =", R_FDR1_Median)) +
  annotate(geom="text", x=-1, y=Rlim -0.1, label=paste("median =", R_FDR1_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Rlim -0.3, label=paste(round(R_FDR1_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Rlim -0.1, label=paste("median =",R_FDR1_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Rlim -0.3, label=paste(round(R_FDR1_perPos, 4) * 100, "%"), hjust=1) +   
  xlim(-1,1) + 
  #ylim(0,3.4) +
  scale_y_continuous(limits = c(0,3.4), breaks=c(0,1,2,3), labels = c("","Roots","Leaves","Flowers")) +
  theme_classic() +
  labs(title = "DE genes FDR < 0.01",
       x = "D index",
       y = "Density"
  ) 

FDR5Density <- cbind(F_FDR5_Density, L_FDR5_Density$y, R_FDR5_Density$y)
names(FDR5Density) <- c("D","Flowers","Leaves","Roots")
head(FDR5Density)

p_fdr5 <- ggplot(data = FDR5Density) +
  geom_line(mapping = aes(x = D, y = Flowers +2),  size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, NA), y = Flowers + 2), fill = "#e55b7e") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Flowers + 2), fill = "#fee2e3") + 
  geom_area(mapping = aes(x = D, y = 2), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Leaves + 1), size=1 )+
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Leaves + 1), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Leaves + 1), fill = "#e6d6c5") + 
  geom_area(mapping = aes(x = D, y = 1), fill ="white") + 
  geom_line(mapping = aes(x = D, y = Roots), size=1 ) +
  geom_area(mapping = aes(x = ifelse(D > 0 , D, 0), y = Roots), fill = "#d5a069") + 
  geom_area(mapping = aes(x = ifelse(D < 0 , D, 0), y = Roots), fill = "#e6d6c5") + 
  geom_segment(aes(x = F_FDR5_Median, y = 2, xend = F_FDR5_Median, yend = 3), color = "gray", linetype = "dashed", size=1) + # Flowers
  geom_segment(aes(x = F_FDR5_Pos, y = 2, xend = F_FDR5_Pos, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = F_FDR5_Neg, y = 2, xend = F_FDR5_Neg, yend = 3), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0.5, y=Flim *1.1, label= "D > 0") +
  annotate(geom="text", x=-0.5, y=Flim *1.1, label= "D < 0") +
  annotate(geom="text", x=0, y=Flim, label=paste("overall median =", F_FDR5_Median)) +
  annotate(geom="text", x=-1, y=Flim, label=paste("median =", F_FDR5_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Flim -0.2, label=paste(round(F_FDR5_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Flim, label=paste("median =",F_FDR5_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Flim -0.2, label=paste(round(F_FDR5_perPos, 4) * 100, "%"), hjust=1) +
  geom_segment(aes(x = L_FDR5_Median, y = 1, xend = L_FDR5_Median, yend = 2), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = L_FDR5_Pos, y = 1, xend = L_FDR5_Pos, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = L_FDR5_Neg, y = 1, xend = L_FDR5_Neg, yend = 2), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Llim -0.1, label=paste("overall median =", L_FDR5_Median)) +
  annotate(geom="text", x=-1, y=Llim -0.1, label=paste("median =", L_FDR5_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Llim -0.3, label=paste(round(L_FDR5_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Llim -0.1, label=paste("median =",L_FDR5_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Llim -0.3, label=paste(round(L_FDR5_perPos, 4) * 100, "%"), hjust=1) +   
  geom_segment(aes(x = R_FDR5_Median, y = 0, xend = R_FDR5_Median, yend = 1), color = "gray", linetype = "dashed", size=1) + # Leaves
  geom_segment(aes(x = R_FDR5_Pos, y = 0, xend = R_FDR5_Pos, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  geom_segment(aes(x = R_FDR5_Neg, y = 0, xend = R_FDR5_Neg, yend = 1), color = "gray34", linetype = "dashed", size=1) +
  annotate(geom="text", x=0, y=Rlim -0.1, label=paste("overall median =", R_FDR5_Median)) +
  annotate(geom="text", x=-1, y=Rlim -0.1, label=paste("median =", R_FDR5_Neg), hjust=0) +
  annotate(geom="text", x=-1, y=Rlim -0.3, label=paste(round(R_FDR5_perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=1, y=Rlim -0.1, label=paste("median =",R_FDR5_Pos), hjust=1) + 
  annotate(geom="text", x=1, y=Rlim -0.3, label=paste(round(R_FDR5_perPos, 4) * 100, "%"), hjust=1) +   
  xlim(-1,1) + 
  #ylim(0,3.4) +
  scale_y_continuous(limits = c(0,3.4), breaks=c(0,1,2,3), labels = c("","Roots","Leaves","Flowers")) +
  theme_classic() +
  labs(title = "DE genes FDR < 0.05",
       x = "D index",
       y = "Density"
  )  

ggarrange(p_fdr1, p_fdr5, p_all, nrow = 1, ncol = 3)  

dim(D_F)[1]
  
ylim <- 1.7
#### Pattern 2
d1 <- ggplot(data = D.Density) +
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Leaves), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = F_All_Median, y = 0, xend = F_All_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = L_All_Median, y = 0, xend = L_All_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = R_All_Median, y = 0, xend = R_All_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=17307)"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=17307)"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=17307)"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "All genes",
       x = "D index (CR)",
       y = "Density"
  ) 

FDR1Density <- cbind(F_FDR1_Density, L_FDR1_Density$y, R_FDR1_Density$y)
names(FDR1Density) <- c("D","Flowers","Leaves","Roots")
head(FDR1Density)


pFDR1 <- ggplot(data = FDR1Density) +
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Roots), fill = "gray87") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Leaves), fill = "gray87") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Flowers), fill = "gray87", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Roots), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Leaves), fill = "white") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#046c46", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#046c46", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "navy",size=1 )+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = ylim * 0.9), color = "gray", linetype = "dashed", size=1) + # Flowers
  #annotate(geom="text", x=0.5, y=ylim, label= "D > 0") +
  #annotate(geom="text", x=-0.5, y=ylim, label= "D < 0") +
  #annotate(geom="text", x=0, y=ylim * 0.9, label=paste("overall median")) +
  #annotate(geom="text", x=-0.7, y=ylim* 0.9, label=paste("median")) +
  #annotate(geom="text", x=0.7, y=ylim* 0.9, label=paste("median")) + 
  #annotate(geom="text", x=-1, y=ylim * 0.8, label=F_FDR1_Neg, color = "red3",hjust=0) +
  #annotate(geom="text", x=-0.7, y=ylim * 0.8, label=L_FDR1_Neg, color = "#046c46") +
  #annotate(geom="text", x=-0.4, y=ylim * 0.8, label=R_FDR1_Neg, color = "#856201",hjust=1) +
  #annotate(geom="text", x=-0.3, y=ylim * 0.8, label=F_FDR1_Median, color = "red3",hjust=0) +
  #annotate(geom="text", x=0, y=ylim * 0.8, label=L_FDR1_Median, color = "#046c46") +
  #annotate(geom="text", x=0.3, y=ylim * 0.8, label=R_FDR1_Median, color = "#856201",hjust=1) +
  #annotate(geom="text", x=0.4, y=ylim * 0.8, label=F_FDR1_Pos, color = "red3",hjust=0) +
  #annotate(geom="text", x=0.7, y=ylim * 0.8, label=L_FDR1_Pos, color = "#046c46") +
  #annotate(geom="text", x=1, y=ylim * 0.8, label=R_FDR1_Pos, color = "#856201",hjust=1) +
  #annotate(geom="text", x=-1, y=ylim * 0.7, label=paste0(round(F_FDR1_perNeg, 4) * 100, "%"), color = "red3", hjust=0) +
  #annotate(geom="text", x=-0.7, y=ylim * 0.7, label=paste0(round(L_FDR1_perNeg, 4) * 100, "%"), color = "#046c46") +
  #annotate(geom="text", x=-0.4, y=ylim * 0.7, label=paste0(round(R_FDR1_perNeg, 4) * 100, "%"), color = "#856201", hjust=1) +
  #annotate(geom="text", x=0.4, y=ylim * 0.7, label=paste0(round(F_FDR1_perPos, 4) * 100, "%"), color = "red3", hjust=0) +
  #annotate(geom="text", x=0.7, y=ylim * 0.7, label=paste0(round(L_FDR1_perPos, 4) * 100, "%"), color = "#046c46") +
  #annotate(geom="text", x=1, y=ylim * 0.7, label=paste0(round(R_FDR1_perPos, 4) * 100, "%"), color = "#856201", hjust=1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "DE genes FDR < 0.01",
       x = "D index",
       y = "Density"
  ) 


FDR5Density <- cbind(F_FDR5_Density, L_FDR5_Density$y, R_FDR5_Density$y)
names(FDR5Density) <- c("D","Flowers","Leaves","Roots")
head(FDR5Density)

d2 <- ggplot(data = FDR5Density) +
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Leaves), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = F_FDR5_Median, y = 0, xend = F_FDR5_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = L_FDR5_Median, y = 0, xend = L_FDR5_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = R_FDR5_Median, y = 0, xend = R_FDR5_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=0.5, y=ylim, label= "D > 0") +
  annotate(geom="text", x=-0.5, y=ylim, label= "D < 0") +
  #annotate(geom="text", x=0.5, y=ylim, label= "D > 0") +
  #annotate(geom="text", x=-0.5, y=ylim, label= "D < 0") +
  #annotate(geom="text", x=0, y=ylim * 0.9, label=paste("overall median")) +
  #annotate(geom="text", x=-0.7, y=ylim* 0.9, label=paste("median")) +
  #annotate(geom="text", x=0.7, y=ylim* 0.9, label=paste("median")) + 
  #annotate(geom="text", x=-1, y=ylim * 0.8, label=F_FDR5_Neg, color = "red3",hjust=0) +
  #annotate(geom="text", x=-0.7, y=ylim * 0.8, label=L_FDR5_Neg, color = "#046c46") +
  #annotate(geom="text", x=-0.4, y=ylim * 0.8, label=R_FDR5_Neg, color = "#856201",hjust=1) +
  #annotate(geom="text", x=-0.3, y=ylim * 0.8, label=F_FDR5_Median, color = "red3",hjust=0) +
  #annotate(geom="text", x=0, y=ylim * 0.8, label=L_FDR5_Median, color = "#046c46") +
  #annotate(geom="text", x=0.3, y=ylim * 0.8, label=R_FDR5_Median, color = "#856201",hjust=1) +
  #annotate(geom="text", x=0.4, y=ylim * 0.8, label=F_FDR5_Pos, color = "red3",hjust=0) +
  #annotate(geom="text", x=0.7, y=ylim * 0.8, label=L_FDR5_Pos, color = "#046c46") +
  #annotate(geom="text", x=1, y=ylim * 0.8, label=R_FDR5_Pos, color = "#856201",hjust=1) +
  #annotate(geom="text", x=-1, y=ylim * 0.7, label=paste0(round(F_FDR5_perNeg, 4) * 100, "%"), color = "red3", hjust=0) +
  #annotate(geom="text", x=-0.7, y=ylim * 0.7, label=paste0(round(L_FDR5_perNeg, 4) * 100, "%"), color = "#046c46") +
  #annotate(geom="text", x=-0.4, y=ylim * 0.7, label=paste0(round(R_FDR5_perNeg, 4) * 100, "%"), color = "#856201", hjust=1) +
  #annotate(geom="text", x=0.4, y=ylim * 0.7, label=paste0(round(F_FDR5_perPos, 4) * 100, "%"), color = "red3", hjust=0) +
  #annotate(geom="text", x=0.7, y=ylim * 0.7, label=paste0(round(L_FDR5_perPos, 4) * 100, "%"), color = "#046c46") +
  #annotate(geom="text", x=1, y=ylim * 0.7, label=paste0(round(R_FDR5_perPos, 4) * 100, "%"), color = "#856201", hjust=1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "DE genes FDR < 0.05",
       x = "D index",
       y = "Density"
  ) 

# CGCR
CGCRDensity <- cbind(F_CGCR_Density, L_CGCR_Density$y, R_CGCR_Density$y)
names(CGCRDensity) <- c("D","Flowers","Leaves","Roots")
head(CGCRDensity)

# find the site that flower overlap with leaf
#      D       Flowers    Leaves    Roots
# -0.29941292 0.5080087 0.5084368 0.4179166

D.overlap <- -0.299

d3 <- ggplot(data = CGCRDensity) +
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D > D.overlap, D, 0), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > D.overlap, D, 0), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D < D.overlap, D, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D < D.overlap, D, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D < D.overlap, D, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D >= D.overlap & D <= 0, D, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = F_CGCR_Median, y = 0, xend = F_CGCR_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = L_CGCR_Median, y = 0, xend = L_CGCR_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = R_CGCR_Median, y = 0, xend = R_CGCR_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCR_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCR_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCR_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCR DE genes",
       x = "D index",
       y = "Density"
  ) 

# CRCO
CRCODensity <- cbind(F_CRCO_Density, L_CRCO_Density$y, R_CRCO_Density$y)
names(CRCODensity) <- c("D","Flowers","Leaves","Roots")
head(CRCODensity)

d4 <- ggplot(data = CRCODensity) +
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, 0), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D < 0, D, 0), y = Leaves), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = F_CRCO_Median, y = 0, xend = F_CRCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = L_CRCO_Median, y = 0, xend = L_CRCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = R_CRCO_Median, y = 0, xend = R_CRCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CRCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CRCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CRCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CRCO DE genes",
       x = "D index",
       y = "Density"
  ) 

# CGCO
CGCODensity <- cbind(F_CGCO_Density, L_CGCO_Density$y, R_CGCO_Density$y)
names(CGCODensity) <- c("D","Flowers","Leaves","Roots")
head(CGCODensity)

# find the site that flower overlap with leaf
#      D       Flowers    Leaves    Roots
# 0.135029354 0.2350747 0.2338749 0.1918482
CGCO.overlap <- 0.135

d5 <- ggplot(data = CGCODensity) +
  geom_area(mapping = aes(x = ifelse(D < 0, D, NA), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D > 0, D, NA), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D < CGCO.overlap, D, NA), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D <= 0, D, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D <= 0, D, NA), y = Leaves), fill = "#f4f4f2") + 
  #geom_area(mapping = aes(x = ifelse(D >= 0 & D <= CGCO.overlap, D, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D >= -0.002 & D <= CGCO.overlap, D, NA), y = Leaves), fill = "white") + 
  geom_line(mapping = aes(x = D, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = F_CGCO_Median, y = 0, xend = F_CGCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = L_CGCO_Median, y = 0, xend = L_CGCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = R_CGCO_Median, y = 0, xend = R_CGCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCO DE genes",
       x = "D index (CR)",
       y = "Density"
  ) 

library(ggpubr)
library("grid")
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(d1, vp = vplayout(1,1))
#print(p2, vp = vplayout(1,2))
print(d3, vp = vplayout(2,1))
print(d4, vp = vplayout(3,1))
print(d5, vp = vplayout(4,1))





