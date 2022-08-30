#  D index pairwise comparison among Cg, Cr, and Co.
## Set work directory and load files
setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
#### This file contains reads counts infromation of each gene in all samples
rawTMM <- read.table("InputData/Capsella_Diploids_ReadsCount_TMM_Filtered.txt", header=TRUE, row.names="Gene")
head(rawTMM)
dim(rawTMM)
##### This file contains information of population and tissue in each sample
Accession <- read.table("InputData/DiploidsPhenotypicFile.txt", header = T)
head(Accession)

### Define populations
CR_F <- c(Accession$species=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$species=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$species=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$species=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$species=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$species=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$species=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$species=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$species=="CO" & Accession$tissue=="root")

geneMis <- 1 
geneKeep <- cbind(  ((rowSums(rawTMM[,CG_F], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CR_F], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CO_F], na.rm = T) >= 1) >= geneMis) &
                      
                      ((rowSums(rawTMM[,CG_L], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CR_L], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CO_L], na.rm = T) >= 1) >= geneMis) &
                      
                      ((rowSums(rawTMM[,CG_R], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CR_R], na.rm = T) >= 1) >= geneMis) &
                      ((rowSums(rawTMM[,CO_R], na.rm = T) >= 1) >= geneMis) 
)
TMM <- rawTMM[geneKeep,]
message(paste("Number of genes before", dim(rawTMM)[1], "\nNumber of genes after", dim(TMM)[1]))

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
getwd()
Random <- read.table("OutputData/Random_Simulation_D.txt", header = T)
head(Random)

density(Random[Random$Value > -0.7 & Random$Value < -0.3,]$Value)


library(ggplot2)
greyScale <- rep(gray.colors(10, start = 0, end = 0.7), 100)     # scale grey color to 10 scales from dark to light, repeat 100 times.

# Plot simulation result
p <- ggplot(Random, aes(x=Value, color = Rep))+
  geom_density() + 
  scale_color_manual(values=greyScale) +   
  xlim(-1,1) + ylim(0,1) +
  theme_classic() +                                 # set theme to classic
  theme(legend.position = "none") +                   # remove legend
  labs( #title = "Fuel economy declines as weight increases",
    #subtitle = "(1973-74)",
    #caption = "Data from the 1974 Motor Trend US magazine.",
    #tag = "D",
    x = " ",
    y = "Density"
    #colour = "Gears"
  )

d <- ggplot_build(p)

head(d$data) # check the density value
Extract <- as.data.frame(d$data) # extract value
head(Extract) 
# Extract$density: density value of y axis
# Extract$group: density groups, from 1 to 1000, 1000random simulations

# Extract values of each simulation within range of -0.7 : -0.3
MyNegative <- Extract[Extract$x > -0.7 & Extract$x < -0.3,]
head(MyNegative)
# calculate the minimum value of each simulation
for(i in 1:1000){
  if (i == 1){
    MyMin <- min(MyNegative[MyNegative$group == i,]$density) # return the minimum density value
    NegativeMin <- MyNegative[((MyNegative$group == i) & (round(MyNegative$density, 7) == round(MyMin,7))),]$x #return the x value, round is important in here
  }else{
    MyMin[i] <- min(MyNegative[MyNegative$group == i,]$density)
    NegativeMin[i] <- MyNegative[((MyNegative$group == i) & (round(MyNegative$density, 7) == round(MyMin[i],7))),]$x #return the x value
  }
}
NegativeMedian <- median(NegativeMin)

# Extract values of each simulation within range of 0.3 : 0.7
MyPositive <- Extract[Extract$x > 0.3 & Extract$x < 0.7,]
head(MyPositive)
# calculate the minimum value of each simulation
for(i in 1:1000){
  if (i == 1){
    MyMin <- min(MyPositive[MyPositive$group == i,]$density) # return the minimum density value
    PositiveMin <- MyPositive[((MyPositive$group == i) & (round(MyPositive$density, 7) == round(MyMin,7))),]$x #return the x value, round is important in here
  }else{
    MyMin[i] <- min(MyPositive[MyPositive$group == i,]$density)
    PositiveMin[i] <- MyPositive[((MyPositive$group == i) & (round(MyPositive$density, 7) == round(MyMin[i],7))),]$x #return the x value
  }
}
PositiveMedian <- median(PositiveMin)

# Calculate percentage of positive, near zero, and negative
head(head(Random)) 
perPos <- dim(Random[Random$Value > PositiveMedian,])[1] /  dim(Random)[1]
perZero <- dim(Random[(Random$Value <= PositiveMedian) & (Random$Value >= NegativeMedian),])[1] / dim(Random)[1]
perNeg <- dim(Random[Random$Value < NegativeMedian,])[1] /  dim(Random)[1]

# plot simulation result

#pdf("~/Desktop/D_10000_Simulations.pdf", width=5, height=4)
pSimulation <- ggplot(Random, aes(x=Value, color = Rep))+
  geom_density() + 
  scale_color_manual(values=greyScale) +   
  xlim(-1,1) + ylim(0,1) +
  annotate(geom="text", x=-1, y=1, label=expression(D[italic(COCR)]~"<"~D[italic(CGCR)]), hjust=0) +
  annotate(geom="text", x=-1, y=1*0.9, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
  annotate(geom="text", x=0, y=1, label=expression(D[italic(COCR)]~"="~D[italic(CGCR)])) +
  annotate(geom="text", x=0, y=1*0.9, label=paste(round(perZero, 4) * 100, "%")) +
  annotate(geom="text", x=1, y=1, label=expression(D[italic(COCR)]~">"~D[italic(CGCR)]), hjust=1) +
  annotate(geom="text", x=1, y=1*0.9, label=paste(round(perPos, 4) * 100, "%"), hjust=1) +
  geom_vline(aes(xintercept=PositiveMedian), color="gray34", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=NegativeMedian), color="gray34", linetype="dashed", size=1) +
  theme_classic() +                                 # set theme to classic
  theme(legend.position = "none") +                   # remove legend
  labs( #title = "Fuel economy declines as weight increases",
    #subtitle = "(1973-74)",
    #caption = "Data from the 1974 Motor Trend US magazine.",
    #tag = "A",
    x = " ",
    y = "Density"
    #colour = "Gears"
  )

#dev.off()

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
  Density <- with(density(D$Value), data.frame(x, y))
  ########################################## 2
  # Seconde loop for All, DE-CGCR/CGCO,CRCO
  for (com in c("All", "CGCR", "CRCO", "CGCO")){
    if (com == "All"){
      N <- dim(D)[1]
      perPos <- dim(D[D$Value > PositiveMedian,])[1] /  dim(D)[1]
      perZero <- dim(D[(D$Value <= PositiveMedian) & (D$Value >= NegativeMedian),])[1] /  dim(D)[1]
      perNeg <- dim(D[D$Value < NegativeMedian,])[1] /  dim(D)[1]
      Mylim <- 1
      de <- "genes" #for plot title
    }else{
      DE <- read.csv(paste("OutputData/DE_", i, "_", com, "_FC2FDR.05.txt", sep = ""), header = T, row.names = 1, sep = "\t")
      keep <- row.names(D) %in% row.names(DE)
      subD <- D[keep,] 
      Density <- with(density(subD$Value), data.frame(x, y))
      N <- dim(subD)[1]
      perPos <- dim(subD[subD$Value > PositiveMedian,])[1] /  dim(subD)[1]
      perZero <- dim(subD[(subD$Value <= PositiveMedian) & (subD$Value >= NegativeMedian),])[1] /  dim(subD)[1]
      perNeg <- dim(subD[subD$Value < NegativeMedian,])[1] /  dim(subD)[1]
      Mylim <- 1.6
      de <- "DE genes"
    } # end if else
    # Plot out
    if (i == "F"){
      p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
        geom_line(size=1 )+
        geom_area(mapping = aes(x = ifelse(x > PositiveMedian , x, 0)), fill = "gray87", size=1) + 
        geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
        geom_area(mapping = aes(x = ifelse(x < NegativeMedian , x, NegativeMedian)), fill = "#ffcecc") + #light red
        geom_vline(aes(xintercept = NegativeMedian), color="gray34", linetype="dashed", size=1) +
        geom_vline(aes(xintercept = PositiveMedian), color="gray34", linetype="dashed", size=1) + 
        xlim(-1,1) + 
        ylim(0,Mylim) +
        annotate(geom="text", x=-1, y=Mylim, label=expression(D[italic(COCR)]~"<"~D[italic(CGCR)]), hjust=0) +
        annotate(geom="text", x=-1, y=Mylim*0.9, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
        annotate(geom="text", x=0, y=Mylim, label=expression(D[italic(COCR)]~"="~D[italic(CGCR)])) +
        annotate(geom="text", x=0, y=Mylim*0.9, label=paste(round(perZero, 4) * 100, "%")) +
        annotate(geom="text", x=1, y=Mylim, label=expression(D[italic(COCR)]~">"~D[italic(CGCR)]), hjust=1) +
        annotate(geom="text", x=1, y=Mylim*0.9, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
        theme_classic() +
        labs(title = paste("Flowers in", com, de,"(n=", N, ")"),
             x = "D index",
             y = "Density"
        ) 
      assign(paste0(i, "_", com, "_p"), p)
    }else{
      if (i == "L"){
        p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
          geom_line(size=1 )+
          geom_area(mapping = aes(x = ifelse(x > PositiveMedian , x, 0)), fill = "gray87", size=1) + 
          geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
          geom_area(mapping = aes(x = ifelse(x < NegativeMedian , x, NegativeMedian)), fill = "#d5f3d3") + #light green
          geom_vline(aes(xintercept = NegativeMedian), color="gray34", linetype="dashed", size=1) +
          geom_vline(aes(xintercept = PositiveMedian), color="gray34", linetype="dashed", size=1) + 
          xlim(-1,1) + ylim(0,Mylim) +
          annotate(geom="text", x=-1, y=Mylim, label=expression(D[italic(COCR)]~"<"~D[italic(CGCR)]), hjust=0) +
          annotate(geom="text", x=-1, y=Mylim * 0.9, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
          annotate(geom="text", x=0, y=Mylim, label=expression(D[italic(COCR)]~"="~D[italic(CGCR)])) +
          annotate(geom="text", x=0, y=Mylim * 0.9, label=paste(round(perZero, 4) * 100, "%")) +
          annotate(geom="text", x=1, y=Mylim, label=expression(D[italic(COCR)]~">"~D[italic(CGCR)]), hjust=1) +
          annotate(geom="text", x=1, y=Mylim * 0.9, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
          theme_classic() +
          labs(title = paste("Leaves in", com, de, "(n=", N, ")"),
               x = "D index",
               y = "Density"
          ) 
        assign(paste0(i, "_", com, "_p"), p)
      }else{
        p <- ggplot(data = Density, mapping = aes(x = x, y = y)) +
          geom_line(size=1 )+
          geom_area(mapping = aes(x = ifelse(x > PositiveMedian , x, 0)), fill = "gray87", size=1) + 
          geom_area(mapping = aes(x = ifelse(x > NegativeMedian & x < PositiveMedian, x, 0)), fill = "white") + 
          geom_area(mapping = aes(x = ifelse(x < NegativeMedian , x, NegativeMedian)), fill = "#f4d8b3") + #light brown
          geom_vline(aes(xintercept = NegativeMedian), color="gray34", linetype="dashed", size=1) +
          geom_vline(aes(xintercept = PositiveMedian), color="gray34", linetype="dashed", size=1) + 
          xlim(-1,1) + ylim(0,Mylim) +
          annotate(geom="text", x=-1, y=Mylim, label=expression(D[italic(COCR)]~"<"~D[italic(CGCR)]), hjust=0) +
          annotate(geom="text", x=-1, y=Mylim*0.9, label=paste(round(perNeg, 4) * 100, "%"), hjust=0) +
          annotate(geom="text", x=0, y=Mylim, label=expression(D[italic(COCR)]~"="~D[italic(CGCR)])) +
          annotate(geom="text", x=0, y=Mylim*0.9, label=paste(round(perZero, 4) * 100, "%")) +
          annotate(geom="text", x=1, y=Mylim, label=expression(D[italic(COCR)]~">"~D[italic(CGCR)]), hjust=1) +
          annotate(geom="text", x=1, y=Mylim*0.9, label=paste(round(perPos, 4) * 100, "%"), hjust=1) + 
          theme_classic() +
          labs(title = paste("Roots in", com, de, "(n=", N, ")"),
               x = "D index",
               y = "Density"
          ) 
        assign(paste0(i, "_", com, "_p"), p)
      }
    } # end if else
  } # end loop 2
} # end loop 1

library(ggpubr)
library("grid")
library(ggthemes) # Load
ggarrange(F_All_p, F_CGCR_p, F_CRCO_p, F_CGCO_p,
          L_All_p, L_CGCR_p, L_CRCO_p, L_CGCO_p,
          R_All_p, R_CGCR_p, R_CRCO_p, R_CGCO_p,
          #labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3, ncol = 4)


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




